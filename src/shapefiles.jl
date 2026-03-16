"""
    list_components(path) -> Vector{NamedTuple}

Print two tables describing the contents of the Shapefile at `path`, then
return a `Vector{NamedTuple}` with one entry per exterior ring (the expanded
table).

**Records table** — one row per DBF record, columns:
`idx | <DBF attributes> | rings`

**Rings table** — one row per exterior ring (the unit produced by
[`read_shapefile`](@ref)), columns:
`idx | ring | <DBF attributes> | n_pts | area | xmin | xmax | ymin | ymax`

- `idx`  — 1-based record index (matches the `AbstractVector{Int}` form of `select`).
- `ring` — 1-based ring index within that record.
- `area` — signed shoelace area in the file's native coordinate units (degrees²
           for geographic CRS; useful for relative size comparison only).
- `xmin/xmax/ymin/ymax` — bounding box of the exterior ring.

The returned `Vector{NamedTuple}` has the same fields as the rings table and can
be used directly with the callable form of `select` in [`read_shapefile`](@ref).
"""
function list_components(path::AbstractString)
  base  = splitext(path)[1]
  table = Shapefile.Table(base * ".shp")
  cols  = _dbf_cols(table)

  # Parse all records.
  record_data = []
  for (i, row) in enumerate(table)
    dbf_vals   = NamedTuple(c => getproperty(row, c) for c in cols)
    shape      = Shapefile.shape(row)
    ring_geoms = isnothing(shape) ? ShapeGeometry[] : _parse_shape(shape)
    push!(record_data, (idx = i, dbf = dbf_vals, rings = ring_geoms))
  end

  isempty(record_data) && return NamedTuple[]

  # --- Records table ---------------------------------------------------------
  n_rec = length(record_data)
  println("── Records ($n_rec) " * "─"^max(1, 50 - ndigits(n_rec) - 12))
  sum_headers  = String["idx"; string.(cols); "rings"]
  sum_str_rows = [[string(rd.idx);
                   [string(rd.dbf[c]) for c in cols];
                   string(length(rd.rings))]
                  for rd in record_data]
  _print_table(sum_headers, sum_str_rows)
  println()

  # --- Rings table -----------------------------------------------------------
  expanded = NamedTuple[]
  for rd in record_data
    for (k, g) in enumerate(rd.rings)
      push!(expanded, _ring_row(rd.idx, k, rd.dbf, g))
    end
  end

  n_rings = length(expanded)
  println("── Rings ($n_rings) " * "─"^max(1, 50 - ndigits(n_rings) - 10))
  if !isempty(expanded)
    geo_cols    = [:n_pts, :area, :xmin, :xmax, :ymin, :ymax]
    exp_headers = String["idx"; "ring"; string.(cols); string.(geo_cols)]
    exp_str_rows = Vector{String}[]
    for r in expanded
      row = String[string(r.idx), string(r.ring)]
      for c in cols; push!(row, string(r[c])); end
      push!(row, string(r.n_pts))
      push!(row, @sprintf("%.4g", r.area))
      push!(row, @sprintf("%.6g", r.xmin)); push!(row, @sprintf("%.6g", r.xmax))
      push!(row, @sprintf("%.6g", r.ymin)); push!(row, @sprintf("%.6g", r.ymax))
      push!(exp_str_rows, row)
    end
    _print_table(exp_headers, exp_str_rows)
  end

  return expanded
end

"""
    read_shapefile(path; select = nothing) -> (Vector{ShapeGeometry}, Union{String,Nothing})

Read a Shapefile and return a vector of geometries plus the raw WKT string
from the `.prj` sidecar file, or `nothing` if no `.prj` is found.

`path` may include or omit the `.shp` extension.

# Keyword arguments
- `select` — restrict which geometries are loaded:
  - `nothing` (default): load all rings from all records.
  - `AbstractVector{Int}`: 1-based **record** indices to keep (as shown in the
    Records table by `list_components`).  All rings within each selected record
    are included.
  - A callable `row -> Bool`: predicate evaluated at **ring level**.  `row` has
    fields `idx`, `ring`, all DBF attribute columns, `n_pts`, `area`, `xmin`,
    `xmax`, `ymin`, `ymax` — matching the Rings table printed by
    `list_components`.  Only rings for which the predicate returns `true` are
    included.  Existing predicates on DBF columns (e.g.
    `row -> row.NAME == "X"`) continue to work unchanged.
- `name_fn` — optional callable `row -> String` evaluated at ring level (same
  `row` as `select`).  When provided, the returned string is stored in the
  `ShapeGeometry.name` field and used as the filename stem when
  `split_components = true`.  When `nothing` (default), names are left empty
  and files are numbered sequentially.

Each exterior ring (plus its holes) becomes one `ShapeGeometry` entry.
"""
function read_shapefile(path::AbstractString; select = nothing, name_fn = nothing)
  base     = splitext(path)[1]
  shp_path = base * ".shp"
  prj_path = base * ".prj"

  source_crs = isfile(prj_path) ? read(prj_path, String) : nothing

  table = Shapefile.Table(shp_path)
  cols  = _dbf_cols(table)
  geoms = ShapeGeometry[]

  for (i, row) in enumerate(table)
    shape = Shapefile.shape(row)
    isnothing(shape) && continue

    if select isa AbstractVector
      i ∈ select || continue
      ring_geoms = _parse_shape(shape)
      if !isnothing(name_fn)
        dbf_vals = NamedTuple(c => getproperty(row, c) for c in cols)
        for (k, g) in enumerate(ring_geoms)
          rrow = _ring_row(i, k, dbf_vals, g)
          push!(geoms, ShapeGeometry(g.exterior, g.holes, name_fn(rrow)))
        end
      else
        append!(geoms, ring_geoms)
      end
    elseif isnothing(select)
      ring_geoms = _parse_shape(shape)
      if !isnothing(name_fn)
        dbf_vals = NamedTuple(c => getproperty(row, c) for c in cols)
        for (k, g) in enumerate(ring_geoms)
          rrow = _ring_row(i, k, dbf_vals, g)
          push!(geoms, ShapeGeometry(g.exterior, g.holes, name_fn(rrow)))
        end
      else
        append!(geoms, ring_geoms)
      end
    else
      # Callable predicate — evaluate at ring level.
      ring_geoms = _parse_shape(shape)
      dbf_vals   = NamedTuple(c => getproperty(row, c) for c in cols)
      for (k, g) in enumerate(ring_geoms)
        rrow = _ring_row(i, k, dbf_vals, g)
        select(rrow) || continue
        name = isnothing(name_fn) ? "" : name_fn(rrow)
        push!(geoms, ShapeGeometry(g.exterior, g.holes, name))
      end
    end
  end

  return geoms, source_crs
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# DBF column names, excluding the pseudo :geometry column.
function _dbf_cols(table)
  filter(!=(:geometry), collect(Tables.columnnames(table)))
end

# Build a ring-level metadata row (used by list_components and read_shapefile).
function _ring_row(idx::Int, ring::Int, dbf_vals::NamedTuple, g::ShapeGeometry)
  pts  = g.exterior.points
  xmin, xmax, ymin, ymax = _bbox(pts)
  merge(
    (idx = idx, ring = ring),
    dbf_vals,
    (n_pts = length(pts),
     area  = _signed_area(pts),
     xmin  = xmin, xmax = xmax,
     ymin  = ymin, ymax = ymax),
  )
end

# Print a table from a header vector and a vector of string-row vectors.
function _print_table(headers::Vector{String}, str_rows::Vector{Vector{String}};
                      max_w::Int = 24)
  widths = [min(max_w, max(length(headers[j]),
                           maximum(length(r[j]) for r in str_rows; init = 0)))
            for j in eachindex(headers)]
  fmt(s, w) = rpad(first(s, w), w)
  hline = join((fmt(headers[j], widths[j]) for j in eachindex(headers)), "  ")
  println(hline)
  println("─"^length(hline))
  for row in str_rows
    println(join((fmt(row[j], widths[j]) for j in eachindex(headers)), "  "))
  end
end

# ---------------------------------------------------------------------------
# Shape parsing (dispatch on Shapefile geometry types)
# ---------------------------------------------------------------------------

function _parse_shape(shape::Shapefile.Polygon)
  return _parse_rings(shape.points, shape.parts)
end

function _parse_shape(shape::Shapefile.Polyline)
  nparts = length(shape.parts)
  npts   = length(shape.points)
  geoms  = ShapeGeometry[]
  for i in 1:nparts
    start_idx = Int(shape.parts[i]) + 1
    end_idx   = i < nparts ? Int(shape.parts[i+1]) : npts
    pts = [(shape.points[j].x, shape.points[j].y) for j in start_idx:end_idx]
    push!(geoms, ShapeGeometry(Contour(pts, false), Contour[]))
  end
  return geoms
end

# Fallback — ignore unsupported geometry types with a warning.
function _parse_shape(shape)
  @warn "Unsupported geometry type: $(typeof(shape)); skipping."
  return ShapeGeometry[]
end

# ---------------------------------------------------------------------------
# Ring parsing and grouping
# ---------------------------------------------------------------------------

# Split raw Shapefile points+parts into classified ShapeGeometry objects.
function _parse_rings(raw_points, parts)
  npts   = length(raw_points)
  nparts = length(parts)

  rings = Vector{NTuple{2,Float64}}[]
  for i in 1:nparts
    start_idx = Int(parts[i]) + 1         # 0-based → 1-based
    end_idx   = i < nparts ? Int(parts[i+1]) : npts
    pts = [(raw_points[j].x, raw_points[j].y) for j in start_idx:end_idx]
    # Shapefiles close rings by repeating the first point — drop the duplicate.
    if length(pts) > 1 && pts[1] == pts[end]
      pop!(pts)
    end
    isempty(pts) && continue
    push!(rings, pts)
  end

  return _group_rings(rings)
end

# Compute the signed area via the shoelace formula.
# Positive → CCW, Negative → CW.
function _signed_area(pts::Vector{NTuple{2,Float64}})
  n = length(pts)
  n < 3 && return 0.0
  area = 0.0
  for i in 1:n
    j = mod1(i + 1, n)
    area += pts[i][1] * pts[j][2] - pts[j][1] * pts[i][2]
  end
  return area / 2
end

# Point-in-polygon test (ray casting).
function _point_in_ring(pt::NTuple{2,Float64}, ring::Vector{NTuple{2,Float64}})
  x, y = pt
  n = length(ring)
  inside = false
  j = n
  for i in 1:n
    xi, yi = ring[i]
    xj, yj = ring[j]
    if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
      inside = !inside
    end
    j = i
  end
  return inside
end

# Group rings into ShapeGeometry objects.
#
# ESRI Shapefile convention (opposite to OGC / math):
#   • CW  rings (negative signed area) → exterior
#   • CCW rings (positive signed area) → hole
#
# We normalise orientation for Gmsh output:
#   • Stored exterior rings → CCW (reversed from ESRI)
#   • Stored hole rings     → CW  (reversed from ESRI)
#
# Each hole is assigned to the smallest exterior ring that contains it.
# A bounding-box pre-filter avoids the expensive point-in-polygon test for
# rings that clearly cannot be contained.
function _group_rings(rings::Vector{Vector{NTuple{2,Float64}}})
  isempty(rings) && return ShapeGeometry[]

  areas    = [_signed_area(r) for r in rings]
  # ESRI: CW (negative) = exterior, CCW (positive) = hole.
  ext_idx  = findall(a -> a <= 0, areas)
  hole_idx = findall(a -> a >  0, areas)

  # Normalise: make exteriors CCW (positive) and holes CW (negative) for Gmsh.
  norm_rings = [copy(r) for r in rings]
  for i in ext_idx
    reverse!(norm_rings[i])   # CW → CCW
  end
  for i in hole_idx
    reverse!(norm_rings[i])   # CCW → CW
  end

  # Pre-compute bounding boxes for all rings.
  bboxes = [_bbox(r) for r in norm_rings]

  # Sort exteriors largest-first (by absolute area) so that iterating in
  # reverse gives the smallest containing exterior first.
  sorted_ext = sort(ext_idx, by = i -> abs(areas[i]), rev = true)

  hole_lists = [Contour[] for _ in sorted_ext]

  for h in hole_idx
    probe    = norm_rings[h][1]
    probe_bb = bboxes[h]
    # Find the smallest exterior that contains the hole.
    for k in length(sorted_ext):-1:1
      e = sorted_ext[k]
      # Quick bbox check before the expensive point-in-polygon test.
      _bbox_contains(bboxes[e], probe_bb) || continue
      if _point_in_ring(probe, norm_rings[e])
        push!(hole_lists[k], Contour(norm_rings[h], true))
        break
      end
    end
    # Holes not contained in any exterior are silently dropped.
  end

  return [ShapeGeometry(Contour(norm_rings[sorted_ext[k]], true), hole_lists[k])
          for k in eachindex(sorted_ext)]
end

# Axis-aligned bounding box as (xmin, xmax, ymin, ymax).
function _bbox(pts::Vector{NTuple{2,Float64}})
  xmin = xmax = pts[1][1]
  ymin = ymax = pts[1][2]
  for (x, y) in pts
    x < xmin && (xmin = x)
    x > xmax && (xmax = x)
    y < ymin && (ymin = y)
    y > ymax && (ymax = y)
  end
  return (xmin, xmax, ymin, ymax)
end

# Return true if bbox `outer` fully contains bbox `inner`.
@inline function _bbox_contains(outer, inner)
  return outer[1] <= inner[1] && outer[2] >= inner[2] &&
         outer[3] <= inner[3] && outer[4] >= inner[4]
end
