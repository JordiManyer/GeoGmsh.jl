"""
    read_shapefile(path) -> (Vector{ShapeGeometry}, Symbol)

Read a Shapefile and return a vector of geometries plus a CRS tag.

`path` may include or omit the `.shp` extension. If a `.prj` sidecar file is
present it is parsed to determine the coordinate reference system, returned as:
- `:degrees`  — geographic CRS (lon/lat in decimal degrees)
- `:meters`   — projected CRS (coordinates already in a linear unit)
- `:unknown`  — no `.prj` found or unrecognised WKT

MultiPolygon records are flattened: each outer ring (plus its holes) becomes a
separate `ShapeGeometry` entry.
"""
function read_shapefile(path::AbstractString)
  base     = splitext(path)[1]
  shp_path = base * ".shp"
  prj_path = base * ".prj"

  crs = isfile(prj_path) ? _read_prj(prj_path) : :unknown

  table = Shapefile.Table(shp_path)
  geoms = ShapeGeometry[]
  for row in table
    shape = Shapefile.shape(row)
    isnothing(shape) && continue
    append!(geoms, _parse_shape(shape))
  end
  return geoms, crs
end

# ---------------------------------------------------------------------------
# CRS detection
# ---------------------------------------------------------------------------

function _read_prj(path::AbstractString) :: Symbol
  wkt = read(path, String)
  if occursin("GEOGCS", wkt) || occursin("Degree", wkt)
    return :degrees
  elseif occursin("PROJCS", wkt)
    return :meters
  else
    return :unknown
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
