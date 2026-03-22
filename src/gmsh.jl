"""
Gmsh geometry file writer and mesh generator.

Two output modes are provided:

- `write_geo`     — writes a human-readable `.geo` script (no Gmsh library
                    call at runtime; open in the GUI or run `gmsh name.geo -2`).
- `generate_mesh` — uses the Gmsh API to build the geometry, mesh it, and
                    write a `.msh` file directly.

Orientation convention (required by Gmsh):
- Exterior curve loops are CCW (positive signed area).
- Hole curve loops are CW (Gmsh subtracts them from the surface).

Because `Geometry2D` already normalises ring orientation (see `shapefiles.jl`),
lines are always written in the stored point order and curve loops always use
positive line indices.
"""

import Gmsh: gmsh

# ============================================================================
# write_geo — text .geo file
# ============================================================================

"""
    write_geo(geoms, name; mesh_size=1.0, mesh_algorithm=nothing,
              split_components=false)

Write `geoms` to a Gmsh `.geo` script.  `name` should be given **without** the
`.geo` extension — it is always appended internally.

When `split_components = false` (default), all geometries are written into a
single file `name.geo`.

When `split_components = true`, one `.geo` file per `Geometry2D` is written
into a directory named `name/`.  Files inside are named `1.geo`, `2.geo`, … .

# Keyword arguments
- `mesh_size`        — characteristic element length assigned to every point.
- `mesh_algorithm`   — optional integer passed as `Mesh.Algorithm` (e.g. 5 =
                       Delaunay, 6 = Frontal-Delaunay, 8 = Frontal-Quad).
                       Omitted from the file when `nothing`.
- `split_components` — write one file per geometry component (default `false`).
"""
function write_geo(
  geoms            :: Vector{Geometry2D},
  name             :: AbstractString;
  mesh_size        :: Real               = 1.0,
  mesh_algorithm   :: Union{Int,Nothing} = nothing,
  split_components :: Bool               = false,
  verbose          :: Bool               = false,
)
  if split_components
    _write_geo_split(geoms, name; mesh_size, mesh_algorithm, verbose)
  else
    _write_geo_single(geoms, name * ".geo"; mesh_size, mesh_algorithm)
    if verbose
      n_pts   = sum(npoints(g.exterior) + sum(npoints(h) for h in g.holes; init=0) for g in geoms)
      n_edges = sum(nedges(g.exterior)  + sum(nedges(h)  for h in g.holes; init=0) for g in geoms)
      println("  Surfaces   : $(length(geoms))")
      println("  Points     : $(_fmt(n_pts))   Edges : $(_fmt(n_edges))")
      println("  Written    : ", name, ".geo")
    end
  end
  return name
end

# ============================================================================
# generate_mesh — .msh file via the Gmsh API
# ============================================================================

"""
    generate_mesh(geoms, name; mesh_size=1.0, mesh_algorithm=nothing,
                  order=1, recombine=false, split_components=false)

Build the geometry with the Gmsh API, generate a 2-D mesh, and write a `.msh`
file.  `name` should be given **without** the `.msh` extension.

When `split_components = true`, one `.msh` file per `Geometry2D` is written
into a directory named `name/`.

# Keyword arguments
- `mesh_size`        — characteristic element length (sets both the min and
                       max bounds passed to Gmsh).
- `mesh_algorithm`   — Gmsh 2-D algorithm tag (e.g. 5 = Delaunay,
                       6 = Frontal-Delaunay, 8 = Frontal-Quad).
                       Uses Gmsh's default when `nothing`.
- `order`            — element order: 1 = linear (default), 2 = quadratic.
- `recombine`        — recombine triangles into quadrilaterals (default `false`).
- `split_components` — write one file per geometry component (default `false`).
"""
function generate_mesh(
  geoms            :: Vector{Geometry2D},
  name             :: AbstractString;
  mesh_size        :: Union{Real,AdaptivityAlgorithm} = 1.0,
  mesh_algorithm   :: Union{Int,Nothing}              = nothing,
  order            :: Int                             = 1,
  recombine        :: Bool                            = false,
  split_components :: Bool                            = false,
  verbose          :: Bool                            = false,
)
  if split_components
    _generate_mesh_split(geoms, name; mesh_size, mesh_algorithm, order, recombine, verbose)
  else
    stats = _generate_mesh_single(geoms, name * ".msh"; mesh_size, mesh_algorithm, order, recombine)
    if verbose
      println("  Nodes      : $(_fmt(stats.nodes))   Elements : $(_fmt(stats.elements))")
      println("  Written    : ", name, ".msh")
    end
  end
  return name
end

# ============================================================================
# Mesh-size helpers (uniform Real vs. AdaptivityAlgorithm)
# ============================================================================

# Characteristic length to pass to Gmsh point entities.
_point_lc(ms::Real)               = Float64(ms)
_point_lc(ms::AdaptivityAlgorithm) = nominal_size(ms)

# Apply mesh-size settings after geometry synchronisation.
# For a plain Real: pin CharacteristicLengthMin == Max == lc.
# For an AdaptivityAlgorithm: set loose bounds and install a background field.
function _setup_mesh_size!(ms::Real, _dem, _curve_tags)
  lc = Float64(ms)
  gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
  gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
end

function _setup_mesh_size!(ms::AdaptivityAlgorithm, dem, curve_tags)
  lc_max = nominal_size(ms)
  gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.0)
  gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_max)
  gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
  gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
  field_tag = apply_adaptivity!(ms, dem, curve_tags)
  isnothing(field_tag) || gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)
end

# ============================================================================
# write_geo internals
# ============================================================================

function _write_geo_single(
  geoms          :: Vector{Geometry2D},
  path           :: AbstractString;
  mesh_size      :: Real,
  mesh_algorithm,
)
  lc = Float64(mesh_size)
  open(path, "w") do io
    _write_header(io, length(geoms), mesh_algorithm)
    pt_id   = 1
    line_id = 1
    loop_id = 1
    surf_id = 1
    for g in geoms
      pt_id, line_id, loop_id, surf_id =
        _write_geometry(io, g, pt_id, line_id, loop_id, surf_id, lc)
    end
  end
end

function _write_geo_split(
  geoms          :: Vector{Geometry2D},
  name           :: AbstractString;
  mesh_size      :: Real,
  mesh_algorithm,
  verbose        :: Bool = false,
)
  mkpath(name)
  n  = length(geoms)
  nd = ndigits(n)
  for (i, g) in enumerate(geoms)
    bname = (isempty(g.name) ? lpad(i, nd, '0') : g.name) * ".geo"
    fname = joinpath(name, bname)
    _write_geo_single([g], fname; mesh_size, mesh_algorithm)
    if verbose
      n_pts   = npoints(g.exterior) + sum(npoints(h) for h in g.holes; init=0)
      n_edges = nedges(g.exterior)  + sum(nedges(h)  for h in g.holes; init=0)
      @printf("  [%*d / %d]  %-*s  %s pts  %s edges\n",
              nd, i, n, nd + 4, bname, _fmt(n_pts), _fmt(n_edges))
    end
  end
  verbose && println("  Written    : ", name, "/")
end

function _write_header(io::IO, ngeoms::Int, mesh_algorithm)
  println(io, "// Generated by ShapefileToGmsh.jl")
  println(io, "// $(ngeoms) surface(s)")
  if !isnothing(mesh_algorithm)
    println(io, "Mesh.Algorithm = $mesh_algorithm;")
  end
  println(io)
end

function _write_geometry(io, g, pt_id, line_id, loop_id, surf_id, lc)
  ext_line_ids, pt_id, line_id =
    _write_contour(io, g.exterior, pt_id, line_id, lc, "exterior")
  ext_loop_id = loop_id
  println(io, "Curve Loop($loop_id) = {$(join(ext_line_ids, ", "))};")
  loop_id += 1

  hole_loop_ids = Int[]
  for (k, hole) in enumerate(g.holes)
    hole_line_ids, pt_id, line_id =
      _write_contour(io, hole, pt_id, line_id, lc, "hole $k")
    cw_ids = reverse(-1 .* hole_line_ids)
    println(io, "Curve Loop($loop_id) = {$(join(cw_ids, ", "))};")
    push!(hole_loop_ids, loop_id)
    loop_id += 1
  end

  all_loops = vcat(ext_loop_id, hole_loop_ids)
  println(io, "Plane Surface($surf_id) = {$(join(all_loops, ", "))};")
  surf_id += 1
  println(io)

  return pt_id, line_id, loop_id, surf_id
end

function _write_contour(
  io      :: IO,
  c       :: Contour,
  pt_id   :: Int,
  line_id :: Int,
  lc      :: Float64,
  label   :: AbstractString,
)
  pts = c.points
  n   = length(pts)
  println(io, "// $label ($n points)")

  first_pt_id = pt_id
  for pt in pts
    @printf(io, "Point(%d) = {%.15g, %.15g, 0, %.15g};\n",
            pt_id, pt[1], pt[2], lc)
    pt_id += 1
  end

  nedges   = c.closed ? n : n - 1
  line_ids = Int[]
  for i in 1:nedges
    p_start = first_pt_id + i - 1
    p_end   = c.closed ? first_pt_id + mod(i, n) : first_pt_id + i
    println(io, "Line($line_id) = {$p_start, $p_end};")
    push!(line_ids, line_id)
    line_id += 1
  end
  println(io)

  return line_ids, pt_id, line_id
end

# ============================================================================
# generate_mesh internals
# ============================================================================

function _generate_mesh_single(
  geoms          :: Vector{Geometry2D},
  path           :: AbstractString;
  mesh_size,
  mesh_algorithm,
  order,
  recombine,
) :: NamedTuple
  lc = _point_lc(mesh_size)
  gmsh.initialize()
  try
    gmsh.option.setNumber("General.Verbosity", 2)   # warnings + errors only
    gmsh.model.add("ShapefileToGmsh")

    pt_id   = 0
    line_id = 0
    loop_id = 0
    surf_id = 0
    for g in geoms
      pt_id, line_id, loop_id, surf_id =
        _add_geometry(g, pt_id, line_id, loop_id, surf_id, lc)
    end

    gmsh.model.geo.synchronize()

    curve_tags = Int[t for (_, t) in gmsh.model.getEntities(1)]
    _setup_mesh_size!(mesh_size, nothing, curve_tags)
    if !isnothing(mesh_algorithm)
      gmsh.option.setNumber("Mesh.Algorithm", Float64(mesh_algorithm))
    end
    if recombine
      gmsh.option.setNumber("Mesh.RecombineAll", 1)
    end

    gmsh.model.mesh.generate(2)
    if order > 1
      gmsh.model.mesh.setOrder(order)
    end

    gmsh.write(path)

    # Collect mesh statistics before finalizing.
    node_tags, _, _   = gmsh.model.mesh.getNodes()
    _, elem_tags, _   = gmsh.model.mesh.getElements(2)   # 2-D elements only
    n_nodes    = length(node_tags)
    n_elements = sum(length(t) for t in elem_tags; init = 0)
    return (; nodes = n_nodes, elements = n_elements)
  finally
    gmsh.finalize()
  end
end

function _generate_mesh_split(
  geoms          :: Vector{Geometry2D},
  name           :: AbstractString;
  mesh_size,
  mesh_algorithm,
  order,
  recombine,
  verbose        :: Bool = false,
)
  mkpath(name)
  n  = length(geoms)
  nd = ndigits(n)
  total_nodes    = 0
  total_elements = 0
  for (i, g) in enumerate(geoms)
    bname = (isempty(g.name) ? lpad(i, nd, '0') : g.name) * ".msh"
    fname = joinpath(name, bname)
    if verbose
      @printf("  [%*d / %d]  %-*s  ", nd, i, n, nd + 4, bname)
      flush(stdout)
    end
    stats = _generate_mesh_single([g], fname; mesh_size, mesh_algorithm, order, recombine)
    total_nodes    += stats.nodes
    total_elements += stats.elements
    if verbose
      @printf("%s nodes  %s elements\n", _fmt(stats.nodes), _fmt(stats.elements))
    end
  end
  if verbose
    println("  Total      : $(_fmt(total_nodes)) nodes  $(_fmt(total_elements)) elements")
    println("  Written    : ", name, "/")
  end
end

# Add one Geometry2D to the active Gmsh model; return updated counters.
function _add_geometry(g, pt_id, line_id, loop_id, surf_id, lc)
  pt_id, line_id, loop_id, ext_loop_id = _add_ring(g.exterior, pt_id, line_id, loop_id, lc)

  hole_loop_ids = Int[]
  for hole in g.holes
    pt_id, line_id, loop_id, hole_loop_id = _add_ring(hole, pt_id, line_id, loop_id, lc)
    push!(hole_loop_ids, hole_loop_id)
  end

  surf_id += 1
  gmsh.model.geo.addPlaneSurface(vcat(ext_loop_id, hole_loop_ids), surf_id)

  return pt_id, line_id, loop_id, surf_id
end

# Add one ring to the active Gmsh model; return updated counters and this loop's tag.
function _add_ring(c, pt_id, line_id, loop_id, lc)
  pts         = c.points
  n           = length(pts)
  first_pt_id = pt_id + 1

  for pt in pts
    pt_id += 1
    gmsh.model.geo.addPoint(pt[1], pt[2], 0.0, lc, pt_id)
  end

  nedges   = c.closed ? n : n - 1
  line_ids = Vector{Int}(undef, nedges)
  for i in 1:nedges
    line_id += 1
    p_start = first_pt_id + i - 1
    p_end   = c.closed ? first_pt_id + mod(i, n) : first_pt_id + i
    gmsh.model.geo.addLine(p_start, p_end, line_id)
    line_ids[i] = line_id
  end

  loop_id += 1
  gmsh.model.geo.addCurveLoop(line_ids, loop_id)

  return pt_id, line_id, loop_id, loop_id
end

# ============================================================================
# 3D overloads — Geometry3D
# ============================================================================

"""
    write_geo(geoms, name; kwargs...)

`Geometry3D` overload: writes boundary points with their terrain elevation
(non-zero z).

When `dem` is provided, interior points are sampled on a regular grid at
spacing `max(mesh_size, DEM pixel size)`, filtered to lie inside each polygon,
and embedded in the surface via `Point{...} In Surface{...};`. This encodes
the full terrain curvature in the `.geo` file so that `gmsh name.geo -2`
produces a terrain-following surface mesh without any additional post-processing.

When `dem` is `nothing` (default), only boundary elevation is written.

Keyword arguments are identical to the 2D version, plus:
- `dem`         — optional [`DEMRaster`](@ref) for interior elevation sampling.
- `nodata_fill` — elevation used for nodata / out-of-bounds cells (default `0.0`).
"""
function write_geo(
  geoms            :: Vector{Geometry3D},
  name             :: AbstractString;
  mesh_size        :: Real                      = 1.0,
  mesh_algorithm   :: Union{Int,Nothing}        = nothing,
  split_components :: Bool                      = false,
  dem              :: Union{DEMRaster,Nothing}  = nothing,
  nodata_fill      :: Real                      = 0.0,
  verbose          :: Bool                      = false,
)
  nf = Float64(nodata_fill)
  if split_components
    mkpath(name)
    n  = length(geoms)
    nd = ndigits(n)
    for (i, g) in enumerate(geoms)
      bname = (isempty(g.base.name) ? lpad(i, nd, '0') : g.base.name) * ".geo"
      _write_geo_single_3d([g], joinpath(name, bname); mesh_size, mesh_algorithm, dem, nodata_fill = nf)
    end
  else
    _write_geo_single_3d(geoms, name * ".geo"; mesh_size, mesh_algorithm, dem, nodata_fill = nf)
    verbose && println("  Written    : ", name, ".geo")
  end
  return name
end

"""
    generate_mesh(geoms, dem, name; kwargs...) -> String

`Geometry3D` overload: generates a terrain-following surface mesh.

Meshes the 2D domain flat, then lifts every node's z-coordinate by sampling
`dem` at its (x, y) position. The DEM and the 2D geometries must be in the
same CRS.

Keyword arguments are identical to the 2D version.
"""
function generate_mesh(
  geoms            :: Vector{Geometry3D},
  dem              :: DEMRaster,
  name             :: AbstractString;
  mesh_size        :: Union{Real,AdaptivityAlgorithm} = 1.0,
  mesh_algorithm   :: Union{Int,Nothing}              = nothing,
  order            :: Int                             = 1,
  recombine        :: Bool                            = false,
  split_components :: Bool                            = false,
  verbose          :: Bool                            = false,
)
  if split_components
    mkpath(name)
    n  = length(geoms)
    nd = ndigits(n)
    for (i, g) in enumerate(geoms)
      bname = (isempty(g.base.name) ? lpad(i, nd, '0') : g.base.name) * ".msh"
      _generate_mesh_single_3d([g], dem, joinpath(name, bname);
        mesh_size, mesh_algorithm, order, recombine)
      verbose && @printf("  [%*d / %d]  %s\n", nd, i, n, bname)
    end
  else
    stats = _generate_mesh_single_3d(geoms, dem, name * ".msh";
      mesh_size, mesh_algorithm, order, recombine)
    if verbose
      println("  Nodes      : $(_fmt(stats.nodes))   Elements : $(_fmt(stats.elements))")
      println("  Written    : ", name, ".msh")
    end
  end
  return name
end

# ---------------------------------------------------------------------------
# 3D internals
# ---------------------------------------------------------------------------

function _write_geo_single_3d(
  geoms          :: Vector{Geometry3D},
  path           :: AbstractString;
  mesh_size      :: Real,
  mesh_algorithm,
  dem            :: Union{DEMRaster,Nothing} = nothing,
  nodata_fill    :: Float64                  = 0.0,
)
  lc = Float64(mesh_size)
  open(path, "w") do io
    _write_header(io, length(geoms), mesh_algorithm)
    pt_id = line_id = loop_id = surf_id = 1
    for g in geoms
      pt_id, line_id, loop_id, surf_id =
        _write_geometry_3d(io, g, pt_id, line_id, loop_id, surf_id, lc, dem, nodata_fill)
    end
  end
end

function _write_geometry_3d(io, g::Geometry3D, pt_id, line_id, loop_id, surf_id, lc,
                             dem::Union{DEMRaster,Nothing}, nodata_fill::Float64)
  ext_line_ids, pt_id, line_id =
    _write_contour_3d(io, g.base.exterior, g.z_exterior, pt_id, line_id, lc, "exterior")
  ext_loop_id = loop_id
  println(io, "Curve Loop($loop_id) = {$(join(ext_line_ids, ", "))};")
  loop_id += 1

  hole_loop_ids = Int[]
  for (k, (hole, z_hole)) in enumerate(zip(g.base.holes, g.z_holes))
    hole_line_ids, pt_id, line_id =
      _write_contour_3d(io, hole, z_hole, pt_id, line_id, lc, "hole $k")
    cw_ids = reverse(-1 .* hole_line_ids)
    println(io, "Curve Loop($loop_id) = {$(join(cw_ids, ", "))};")
    push!(hole_loop_ids, loop_id)
    loop_id += 1
  end

  all_loops = vcat(ext_loop_id, hole_loop_ids)
  # Use Plane Surface; Gmsh will project to best-fit plane for non-flat domains.
  println(io, "Plane Surface($surf_id) = {$(join(all_loops, ", "))};")

  if !isnothing(dem)
    interior_pt_ids, pt_id = _write_interior_dem_points!(
      io, g.base, dem, pt_id, lc, nodata_fill)
    if !isempty(interior_pt_ids)
      println(io, "Point{$(join(interior_pt_ids, ", "))} In Surface{$surf_id};")
    end
  end

  surf_id += 1
  println(io)

  return pt_id, line_id, loop_id, surf_id
end

# Sample DEM on a regular grid inside the polygon, write Point entries, return ids.
function _write_interior_dem_points!(
  io          :: IO,
  g           :: Geometry2D,
  dem         :: DEMRaster,
  pt_id_start :: Int,
  lc          :: Float64,
  nodata_fill :: Float64,
) :: Tuple{Vector{Int}, Int}
  # Grid spacing: no finer than DEM pixel size
  spacing = max(lc, abs(dem.transform[2]))

  xs = [p[1] for p in g.exterior.points]
  ys = [p[2] for p in g.exterior.points]
  x_min, x_max = minimum(xs), maximum(xs)
  y_min, y_max = minimum(ys), maximum(ys)

  # Collect candidate grid points strictly inside the polygon
  interior_pts = NTuple{2,Float64}[]
  x = x_min + spacing
  while x < x_max
    y = y_min + spacing
    while y < y_max
      if _point_in_polygon((x, y), g.exterior, g.holes)
        push!(interior_pts, (x, y))
      end
      y += spacing
    end
    x += spacing
  end

  isempty(interior_pts) && return Int[], pt_id_start

  z_vals = sample_elevation(interior_pts, dem; nodata_fill)

  println(io, "// interior DEM points ($(length(interior_pts)))")
  pt_ids = Int[]
  pt_id  = pt_id_start
  for ((xi, yi), zi) in zip(interior_pts, z_vals)
    @printf(io, "Point(%d) = {%.15g, %.15g, %.15g, %.15g};\n", pt_id, xi, yi, zi, lc)
    push!(pt_ids, pt_id)
    pt_id += 1
  end
  println(io)

  return pt_ids, pt_id
end

function _point_in_polygon(pt::NTuple{2,Float64}, exterior::Contour, holes::Vector{Contour})
  _ray_cast(pt, exterior.points) || return false
  for hole in holes
    _ray_cast(pt, hole.points) && return false
  end
  return true
end

function _ray_cast(pt::NTuple{2,Float64}, ring::Vector{NTuple{2,Float64}})
  x, y   = pt
  n      = length(ring)
  inside = false
  j      = n
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

function _write_contour_3d(
  io      :: IO,
  c       :: Contour,
  z_vals  :: Vector{Float64},
  pt_id   :: Int,
  line_id :: Int,
  lc      :: Float64,
  label   :: AbstractString,
)
  pts = c.points
  n   = length(pts)
  println(io, "// $label ($n points, terrain z)")

  first_pt_id = pt_id
  for (pt, z) in zip(pts, z_vals)
    @printf(io, "Point(%d) = {%.15g, %.15g, %.15g, %.15g};\n",
            pt_id, pt[1], pt[2], z, lc)
    pt_id += 1
  end

  nedges_n = c.closed ? n : n - 1
  line_ids = Int[]
  for i in 1:nedges_n
    p_start = first_pt_id + i - 1
    p_end   = c.closed ? first_pt_id + mod(i, n) : first_pt_id + i
    println(io, "Line($line_id) = {$p_start, $p_end};")
    push!(line_ids, line_id)
    line_id += 1
  end
  println(io)

  return line_ids, pt_id, line_id
end

# ============================================================================
# generate_mesh_volume — volumetric (tetrahedral) mesh via prism extrusion
# ============================================================================

"""
    generate_mesh_volume(geoms, dem, name; depth=nothing, z_bottom=nothing,
                         z_top=nothing, mesh_size=1.0, mesh_algorithm=nothing,
                         split_components=false, verbose=false)

Generate a volumetric tetrahedral mesh by extruding the terrain surface to
flat bounding planes.

The extrusion direction is controlled by `depth`:
- `depth = d` (plain `Real`) — extrude **downward** by `d` below `min(z_terrain)`.
- `depth = (d_below, d_above)` (2-tuple) — extrude **both** directions
  simultaneously.  The terrain surface becomes an interior **Interface** physical
  group shared by **Volume_Below** and **Volume_Above**.
- `depth = (0.0, h)` — extrude **upward** only.

`z_bottom` and `z_top` pin the flat planes to absolute elevations and override
the corresponding component of `depth`.

Physical groups written to the `.msh` file:

| Scenario | Groups |
|---|---|
| Downward only | `Volume`, `Top` (terrain), `Bottom` (flat), `Sides` |
| Upward only   | `Volume`, `Top` (flat), `Bottom` (terrain), `Sides` |
| Both          | `Volume_Below`, `Volume_Above`, `Interface` (terrain), `Top`, `Bottom`, `Sides` |

# Arguments
- `geoms` — `Vector{Geometry3D}` (only the 2D base geometry is used).
- `dem`   — [`DEMRaster`](@ref) for elevation sampling.
- `name`  — output path **without** the `.msh` extension.

# Keyword arguments
- `depth`            — see above.
- `z_bottom`         — absolute bottom elevation; overrides the below component.
- `z_top`            — absolute top elevation; overrides the above component.
- `mesh_size`        — characteristic element length (default `1.0`).
- `mesh_algorithm`   — Gmsh 2D algorithm tag.
- `split_components` — write one `.msh` per geometry component (default `false`).
- `verbose`          — print progress (default `false`).
"""
function generate_mesh_volume(
  geoms            :: Vector{Geometry3D},
  dem              :: DEMRaster,
  name             :: AbstractString;
  depth            :: Union{Real,Tuple{Real,Real},Nothing} = nothing,
  z_bottom         :: Union{Real,Nothing}                  = nothing,
  z_top            :: Union{Real,Nothing}                  = nothing,
  mesh_size        :: Union{Real,AdaptivityAlgorithm}      = 1.0,
  mesh_algorithm   :: Union{Int,Nothing}                   = nothing,
  split_components :: Bool                                 = false,
  verbose          :: Bool                                 = false,
)
  # Parse depth → (d_below, d_above)
  d_below, d_above = if isnothing(depth)
    0.0, 0.0
  elseif depth isa Real
    Float64(depth), 0.0
  else
    Float64(depth[1]), Float64(depth[2])
  end

  has_below = d_below > 0.0 || !isnothing(z_bottom)
  has_above = d_above > 0.0 || !isnothing(z_top)
  has_below || has_above ||
    error("Specify at least one of: `depth`, `z_bottom`, or `z_top`.")

  z_bot_f = isnothing(z_bottom) ? nothing : Float64(z_bottom)
  z_top_f = isnothing(z_top)    ? nothing : Float64(z_top)

  _mesh_vol(gs, path) = _generate_mesh_volume_single(gs, dem, path;
    mesh_size, mesh_algorithm,
    d_below, d_above, z_bottom = z_bot_f, z_top = z_top_f)

  if split_components
    mkpath(name)
    n  = length(geoms)
    nd = ndigits(n)
    total_nodes = total_tets = 0
    for (i, g) in enumerate(geoms)
      bname = (isempty(g.base.name) ? lpad(i, nd, '0') : g.base.name) * ".msh"
      fpath = joinpath(name, bname)
      verbose && @printf("  [%*d / %d]  %-*s  ", nd, i, n, nd + 4, bname)
      stats = _mesh_vol([g], fpath)
      total_nodes += stats.nodes
      total_tets  += stats.elements
      verbose && @printf("%s nodes  %s tets\n", _fmt(stats.nodes), _fmt(stats.elements))
    end
    if verbose
      println("  Total      : $(_fmt(total_nodes)) nodes  $(_fmt(total_tets)) tets")
      println("  Written    : ", name, "/")
    end
  else
    stats = _mesh_vol(geoms, name * ".msh")
    if verbose
      println("  Nodes      : $(_fmt(stats.nodes))   Tets : $(_fmt(stats.elements))")
      println("  Written    : ", name, ".msh")
    end
  end
  return name
end

# ── Shared step 1–3: 2D mesh, DEM sampling, boundary edges ─────────────────

function _generate_mesh_volume_single(
  geoms        :: Vector{Geometry3D},
  dem          :: DEMRaster,
  path         :: AbstractString;
  mesh_size,
  mesh_algorithm,
  d_below      :: Float64,
  d_above      :: Float64,
  z_bottom     :: Union{Float64,Nothing},
  z_top        :: Union{Float64,Nothing},
) :: NamedTuple
  lc = _point_lc(mesh_size)

  # ── Step 1: flat 2D triangle mesh ──────────────────────────────────────────
  node_tags_2d   = Vector{Int}()
  coords_2d      = Vector{Float64}()
  tri_nodes_flat = Vector{Int}()

  gmsh.initialize()
  try
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("_2d_tmp")
    pt_id = line_id = loop_id = surf_id = 0
    for g in geoms
      pt_id, line_id, loop_id, surf_id =
        _add_geometry(g.base, pt_id, line_id, loop_id, surf_id, lc)
    end
    gmsh.model.geo.synchronize()
    curve_tags = Int[t for (_, t) in gmsh.model.getEntities(1)]
    _setup_mesh_size!(mesh_size, dem, curve_tags)
    !isnothing(mesh_algorithm) &&
      gmsh.option.setNumber("Mesh.Algorithm", Float64(mesh_algorithm))
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    append!(node_tags_2d, tags)
    append!(coords_2d,    coords)
    etypes, _, enodes = gmsh.model.mesh.getElements(2)
    tri_idx = findfirst(==(2), etypes)
    !isnothing(tri_idx) && append!(tri_nodes_flat, enodes[tri_idx])
  finally
    gmsh.finalize()
  end

  n_nodes    = length(node_tags_2d)
  n_tris     = length(tri_nodes_flat) ÷ 3
  tag_to_idx = Dict{Int,Int}(tag => i for (i, tag) in enumerate(node_tags_2d))

  # ── Step 2: DEM sampling + elevation resolution ────────────────────────────
  pts_2d = NTuple{2,Float64}[(coords_2d[3i-2], coords_2d[3i-1]) for i in 1:n_nodes]
  z_surf = sample_elevation(pts_2d, dem)

  has_below = d_below > 0.0 || !isnothing(z_bottom)
  has_above = d_above > 0.0 || !isnothing(z_top)
  z_bot_abs = has_below ? (isnothing(z_bottom) ? minimum(z_surf) - d_below : z_bottom) : nothing
  z_top_abs = has_above ? (isnothing(z_top)    ? maximum(z_surf) + d_above : z_top)    : nothing

  # ── Step 3: directed boundary edges (for side-wall triangles) ──────────────
  directed_edges = Set{NTuple{2,Int}}()
  for k in 1:n_tris
    i0 = tag_to_idx[tri_nodes_flat[3k-2]]
    i1 = tag_to_idx[tri_nodes_flat[3k-1]]
    i2 = tag_to_idx[tri_nodes_flat[3k  ]]
    push!(directed_edges, (i0, i1), (i1, i2), (i2, i0))
  end
  boundary_edges = NTuple{2,Int}[(a, b) for (a, b) in directed_edges if (b, a) ∉ directed_edges]

  # ── Step 4: assemble ───────────────────────────────────────────────────────
  if has_below && has_above
    return _volume_bidirectional(path, coords_2d, z_surf, z_bot_abs, z_top_abs,
                                 tag_to_idx, tri_nodes_flat, n_nodes, n_tris,
                                 boundary_edges)
  else
    z_flat    = has_below ? z_bot_abs : z_top_abs
    is_upward = !has_below
    return _volume_unidirectional(path, coords_2d, z_surf, z_flat, is_upward,
                                  tag_to_idx, tri_nodes_flat, n_nodes, n_tris,
                                  boundary_edges)
  end
end

# ── Unidirectional extrusion (down OR up) ────────────────────────────────────
#
# Node layout:  surf (terrain) = 1..n_nodes,  flat = n_nodes+1..2*n_nodes
#
# DOWN (is_upward=false): flat is below terrain
#   • Tets:  b=flat, t=surf
#   • Top    = terrain (original winding, normal ↑)
#   • Bottom = flat    (reversed winding,  normal ↓)
#   • Sides: T1=(surf_a, flat_a, flat_b)  T2=(surf_a, flat_b, surf_b)
#
# UP (is_upward=true): flat is above terrain
#   • Tets:  b=surf, t=flat
#   • Top    = flat    (original winding, normal ↑)
#   • Bottom = terrain (reversed winding, normal ↓)
#   • Sides: T1=(surf_a, surf_b, flat_b)  T2=(surf_a, flat_b, flat_a)

function _volume_unidirectional(
  path, coords_2d, z_surf, z_flat, is_upward,
  tag_to_idx, tri_nodes_flat, n_nodes, n_tris, boundary_edges,
) :: NamedTuple
  n_tets       = 3 * n_tris
  n_bnd        = length(boundary_edges)
  n_side_tris  = 2 * n_bnd

  gmsh.initialize()
  try
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("GeoGmsh3DVolume")

    vol_tag      = gmsh.model.addDiscreteEntity(3)
    surf_ent     = gmsh.model.addDiscreteEntity(2)   # terrain surface entity
    flat_ent     = gmsh.model.addDiscreteEntity(2)   # flat surface entity
    side_ent     = gmsh.model.addDiscreteEntity(2)

    # All nodes to the volume entity
    all_node_tags = collect(1:2*n_nodes)
    all_coords    = Vector{Float64}(undef, 6*n_nodes)
    for i in 1:n_nodes
      x = coords_2d[3i-2];  y = coords_2d[3i-1];  fi = n_nodes + i
      all_coords[3i-2]  = x;  all_coords[3i-1]  = y;  all_coords[3i  ] = z_surf[i]
      all_coords[3fi-2] = x;  all_coords[3fi-1] = y;  all_coords[3fi ] = z_flat
    end
    gmsh.model.mesh.addNodes(3, vol_tag, all_node_tags, all_coords)

    # Tetrahedra — prism split: 3 tets per triangle
    #   DOWN: b = flat (n_nodes+i),  t = surf (i)
    #   UP:   b = surf (i),          t = flat (n_nodes+i)
    tet_tags  = collect(1:n_tets)
    tet_nodes = Vector{Int}(undef, 4*n_tets)
    for k in 1:n_tris
      si0 = tag_to_idx[tri_nodes_flat[3k-2]]
      si1 = tag_to_idx[tri_nodes_flat[3k-1]]
      si2 = tag_to_idx[tri_nodes_flat[3k  ]]
      b0, b1, b2 = is_upward ? (si0, si1, si2) : (n_nodes+si0, n_nodes+si1, n_nodes+si2)
      t0, t1, t2 = is_upward ? (n_nodes+si0, n_nodes+si1, n_nodes+si2) : (si0, si1, si2)
      base = 12*(k-1)
      tet_nodes[base+ 1]=b0;  tet_nodes[base+ 2]=b1
      tet_nodes[base+ 3]=b2;  tet_nodes[base+ 4]=t0
      tet_nodes[base+ 5]=t0;  tet_nodes[base+ 6]=b1
      tet_nodes[base+ 7]=b2;  tet_nodes[base+ 8]=t2
      tet_nodes[base+ 9]=t0;  tet_nodes[base+10]=t1
      tet_nodes[base+11]=b1;  tet_nodes[base+12]=t2
    end
    gmsh.model.mesh.addElements(3, vol_tag, [4], [tet_tags], [tet_nodes])

    # Surface triangles
    #   "Top" = higher-z surface  → original winding (normal ↑)
    #   "Bottom" = lower-z surface → reversed winding (normal ↓)
    surf_elem_tags = collect(n_tets+1 : n_tets+n_tris)
    flat_elem_tags = collect(n_tets+n_tris+1 : n_tets+2*n_tris)
    surf_nodes_v   = Vector{Int}(undef, 3*n_tris)   # original winding
    flat_nodes_v   = Vector{Int}(undef, 3*n_tris)   # original winding before reversal
    for k in 1:n_tris
      si0 = tag_to_idx[tri_nodes_flat[3k-2]]
      si1 = tag_to_idx[tri_nodes_flat[3k-1]]
      si2 = tag_to_idx[tri_nodes_flat[3k  ]]
      surf_nodes_v[3k-2]=si0;         surf_nodes_v[3k-1]=si1;         surf_nodes_v[3k  ]=si2
      flat_nodes_v[3k-2]=n_nodes+si0; flat_nodes_v[3k-1]=n_nodes+si1; flat_nodes_v[3k  ]=n_nodes+si2
    end
    # DOWN: surf=Top (orig), flat=Bottom (reversed)
    # UP:   flat=Top (orig), surf=Bottom (reversed)
    if is_upward
      # Reverse surf winding for "Bottom"
      for k in 1:n_tris
        surf_nodes_v[3k-1], surf_nodes_v[3k] = surf_nodes_v[3k], surf_nodes_v[3k-1]
      end
      top_ent, top_tags_v, top_nodes_arr = flat_ent, flat_elem_tags, flat_nodes_v
      bot_ent, bot_tags_v, bot_nodes_arr = surf_ent, surf_elem_tags, surf_nodes_v
      top_phys, bot_phys = "Top", "Bottom"
    else
      # Reverse flat winding for "Bottom"
      for k in 1:n_tris
        flat_nodes_v[3k-1], flat_nodes_v[3k] = flat_nodes_v[3k], flat_nodes_v[3k-1]
      end
      top_ent, top_tags_v, top_nodes_arr = surf_ent, surf_elem_tags, surf_nodes_v
      bot_ent, bot_tags_v, bot_nodes_arr = flat_ent, flat_elem_tags, flat_nodes_v
      top_phys, bot_phys = "Top", "Bottom"
    end
    gmsh.model.mesh.addElements(2, top_ent, [2], [top_tags_v], [top_nodes_arr])
    gmsh.model.mesh.addElements(2, bot_ent, [2], [bot_tags_v], [bot_nodes_arr])

    # Side-wall triangles (2 per boundary edge, outward normal)
    side_elem_tags = collect(n_tets+2*n_tris+1 : n_tets+2*n_tris+n_side_tris)
    side_nodes_v   = Vector{Int}(undef, 3*n_side_tris)
    for (j, (a, b)) in enumerate(boundary_edges)
      fa = n_nodes+a;  fb = n_nodes+b
      base = 6*(j-1)
      if is_upward
        # surf=lower, flat=upper  →  T1:(a,b,fb)  T2:(a,fb,fa)
        side_nodes_v[base+1]=a;  side_nodes_v[base+2]=b;  side_nodes_v[base+3]=fb
        side_nodes_v[base+4]=a;  side_nodes_v[base+5]=fb; side_nodes_v[base+6]=fa
      else
        # surf=upper, flat=lower  →  T1:(a,fa,fb)  T2:(a,fb,b)
        side_nodes_v[base+1]=a;  side_nodes_v[base+2]=fa; side_nodes_v[base+3]=fb
        side_nodes_v[base+4]=a;  side_nodes_v[base+5]=fb; side_nodes_v[base+6]=b
      end
    end
    gmsh.model.mesh.addElements(2, side_ent, [2], [side_elem_tags], [side_nodes_v])

    gmsh.model.addPhysicalGroup(3, [vol_tag],  1, "Volume")
    gmsh.model.addPhysicalGroup(2, [top_ent],  2, top_phys)
    gmsh.model.addPhysicalGroup(2, [bot_ent],  3, bot_phys)
    gmsh.model.addPhysicalGroup(2, [side_ent], 4, "Sides")

    gmsh.write(path)
    return (; nodes = 2*n_nodes, elements = n_tets)
  finally
    gmsh.finalize()
  end
end

# ── Bidirectional extrusion (down AND up) ────────────────────────────────────
#
# Node layout:
#   surf (terrain) = 1..n_nodes        z = z_surf[i]
#   bot flat       = n_nodes+1..2n     z = z_bot  (constant)
#   top flat       = 2*n_nodes+1..3n   z = z_top  (constant)
#
# Physical groups: Volume_Below, Volume_Above, Interface (terrain),
#                  Top (flat top), Bottom (flat bot), Sides

function _volume_bidirectional(
  path, coords_2d, z_surf, z_bot, z_top,
  tag_to_idx, tri_nodes_flat, n_nodes, n_tris, boundary_edges,
) :: NamedTuple
  n_tets_half  = 3 * n_tris       # tets in each half-volume
  n_tets       = 2 * n_tets_half
  n_bnd        = length(boundary_edges)
  n_side_each  = 2 * n_bnd        # side triangles per half-volume

  gmsh.initialize()
  try
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("GeoGmsh3DVolumeBI")

    vol_below = gmsh.model.addDiscreteEntity(3)
    vol_above = gmsh.model.addDiscreteEntity(3)
    iface_ent = gmsh.model.addDiscreteEntity(2)   # Interface (terrain)
    bot_ent   = gmsh.model.addDiscreteEntity(2)   # flat bottom
    top_ent   = gmsh.model.addDiscreteEntity(2)   # flat top
    side_ent  = gmsh.model.addDiscreteEntity(2)   # all side walls

    # All nodes stored on vol_below; referenced by all entities
    all_node_tags = collect(1:3*n_nodes)
    all_coords    = Vector{Float64}(undef, 9*n_nodes)
    for i in 1:n_nodes
      x = coords_2d[3i-2];  y = coords_2d[3i-1]
      bi = n_nodes   + i;    ti = 2*n_nodes + i
      all_coords[3i-2]  = x;  all_coords[3i-1]  = y;  all_coords[3i  ] = z_surf[i]
      all_coords[3bi-2] = x;  all_coords[3bi-1] = y;  all_coords[3bi ] = z_bot
      all_coords[3ti-2] = x;  all_coords[3ti-1] = y;  all_coords[3ti ] = z_top
    end
    gmsh.model.mesh.addNodes(3, vol_below, all_node_tags, all_coords)

    # Tets below: b=bot_flat, t=surf(terrain)
    tet_below_tags  = collect(1:n_tets_half)
    tet_below_nodes = Vector{Int}(undef, 4*n_tets_half)
    for k in 1:n_tris
      si0=tag_to_idx[tri_nodes_flat[3k-2]]
      si1=tag_to_idx[tri_nodes_flat[3k-1]]
      si2=tag_to_idx[tri_nodes_flat[3k  ]]
      b0=n_nodes+si0; b1=n_nodes+si1; b2=n_nodes+si2
      t0=si0;         t1=si1;         t2=si2
      base=12*(k-1)
      tet_below_nodes[base+ 1]=b0;  tet_below_nodes[base+ 2]=b1
      tet_below_nodes[base+ 3]=b2;  tet_below_nodes[base+ 4]=t0
      tet_below_nodes[base+ 5]=t0;  tet_below_nodes[base+ 6]=b1
      tet_below_nodes[base+ 7]=b2;  tet_below_nodes[base+ 8]=t2
      tet_below_nodes[base+ 9]=t0;  tet_below_nodes[base+10]=t1
      tet_below_nodes[base+11]=b1;  tet_below_nodes[base+12]=t2
    end
    gmsh.model.mesh.addElements(3, vol_below, [4], [tet_below_tags], [tet_below_nodes])

    # Tets above: b=surf(terrain), t=top_flat
    tet_above_tags  = collect(n_tets_half+1:n_tets)
    tet_above_nodes = Vector{Int}(undef, 4*n_tets_half)
    for k in 1:n_tris
      si0=tag_to_idx[tri_nodes_flat[3k-2]]
      si1=tag_to_idx[tri_nodes_flat[3k-1]]
      si2=tag_to_idx[tri_nodes_flat[3k  ]]
      b0=si0;              b1=si1;              b2=si2
      t0=2*n_nodes+si0;    t1=2*n_nodes+si1;    t2=2*n_nodes+si2
      base=12*(k-1)
      tet_above_nodes[base+ 1]=b0;  tet_above_nodes[base+ 2]=b1
      tet_above_nodes[base+ 3]=b2;  tet_above_nodes[base+ 4]=t0
      tet_above_nodes[base+ 5]=t0;  tet_above_nodes[base+ 6]=b1
      tet_above_nodes[base+ 7]=b2;  tet_above_nodes[base+ 8]=t2
      tet_above_nodes[base+ 9]=t0;  tet_above_nodes[base+10]=t1
      tet_above_nodes[base+11]=b1;  tet_above_nodes[base+12]=t2
    end
    gmsh.model.mesh.addElements(3, vol_above, [4], [tet_above_tags], [tet_above_nodes])

    # Interface: terrain surface, original winding (normal ↑, faces Volume_Below outward)
    iface_tags = collect(n_tets+1:n_tets+n_tris)
    iface_nodes = Vector{Int}(undef, 3*n_tris)
    for k in 1:n_tris
      iface_nodes[3k-2]=tag_to_idx[tri_nodes_flat[3k-2]]
      iface_nodes[3k-1]=tag_to_idx[tri_nodes_flat[3k-1]]
      iface_nodes[3k  ]=tag_to_idx[tri_nodes_flat[3k  ]]
    end
    gmsh.model.mesh.addElements(2, iface_ent, [2], [iface_tags], [iface_nodes])

    # Bottom: flat bottom, reversed winding (normal ↓)
    bot_elem_tags = collect(n_tets+n_tris+1:n_tets+2*n_tris)
    bot_nodes_v   = Vector{Int}(undef, 3*n_tris)
    for k in 1:n_tris
      si0=tag_to_idx[tri_nodes_flat[3k-2]]
      si1=tag_to_idx[tri_nodes_flat[3k-1]]
      si2=tag_to_idx[tri_nodes_flat[3k  ]]
      bot_nodes_v[3k-2]=n_nodes+si0; bot_nodes_v[3k-1]=n_nodes+si2; bot_nodes_v[3k]=n_nodes+si1
    end
    gmsh.model.mesh.addElements(2, bot_ent, [2], [bot_elem_tags], [bot_nodes_v])

    # Top: flat top, original winding (normal ↑)
    top_elem_tags = collect(n_tets+2*n_tris+1:n_tets+3*n_tris)
    top_nodes_v   = Vector{Int}(undef, 3*n_tris)
    for k in 1:n_tris
      si0=tag_to_idx[tri_nodes_flat[3k-2]]
      si1=tag_to_idx[tri_nodes_flat[3k-1]]
      si2=tag_to_idx[tri_nodes_flat[3k  ]]
      top_nodes_v[3k-2]=2*n_nodes+si0; top_nodes_v[3k-1]=2*n_nodes+si1; top_nodes_v[3k]=2*n_nodes+si2
    end
    gmsh.model.mesh.addElements(2, top_ent, [2], [top_elem_tags], [top_nodes_v])

    # Sides: lower half (surf→bot) + upper half (surf→top)
    n_side_total  = 2 * n_side_each
    side_elem_tags = collect(n_tets+3*n_tris+1 : n_tets+3*n_tris+n_side_total)
    side_nodes_v   = Vector{Int}(undef, 3*n_side_total)
    for (j, (a, b)) in enumerate(boundary_edges)
      ba=n_nodes+a;      bb=n_nodes+b
      ta=2*n_nodes+a;    tb=2*n_nodes+b
      # Lower (surf→bot): T1=(a,ba,bb)  T2=(a,bb,b)
      lo = 6*(j-1)
      side_nodes_v[lo+1]=a;  side_nodes_v[lo+2]=ba; side_nodes_v[lo+3]=bb
      side_nodes_v[lo+4]=a;  side_nodes_v[lo+5]=bb; side_nodes_v[lo+6]=b
      # Upper (surf→top): T1=(a,b,tb)  T2=(a,tb,ta)
      hi = 3*n_side_each + 6*(j-1)
      side_nodes_v[hi+1]=a;  side_nodes_v[hi+2]=b;  side_nodes_v[hi+3]=tb
      side_nodes_v[hi+4]=a;  side_nodes_v[hi+5]=tb; side_nodes_v[hi+6]=ta
    end
    gmsh.model.mesh.addElements(2, side_ent, [2], [side_elem_tags], [side_nodes_v])

    gmsh.model.addPhysicalGroup(3, [vol_below], 1, "Volume_Below")
    gmsh.model.addPhysicalGroup(3, [vol_above], 2, "Volume_Above")
    gmsh.model.addPhysicalGroup(2, [iface_ent], 3, "Interface")
    gmsh.model.addPhysicalGroup(2, [bot_ent],   4, "Bottom")
    gmsh.model.addPhysicalGroup(2, [top_ent],   5, "Top")
    gmsh.model.addPhysicalGroup(2, [side_ent],  6, "Sides")

    gmsh.write(path)
    return (; nodes = 3*n_nodes, elements = n_tets)
  finally
    gmsh.finalize()
  end
end

function _generate_mesh_single_3d(
  geoms          :: Vector{Geometry3D},
  dem            :: DEMRaster,
  path           :: AbstractString;
  mesh_size,
  mesh_algorithm,
  order,
  recombine,
) :: NamedTuple
  lc = _point_lc(mesh_size)
  gmsh.initialize()
  try
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("GeoGmsh3D")

    # Build flat 2D geometry (z=0) from the 2D base of each Geometry3D.
    pt_id = line_id = loop_id = surf_id = 0
    for g in geoms
      pt_id, line_id, loop_id, surf_id =
        _add_geometry(g.base, pt_id, line_id, loop_id, surf_id, lc)
    end

    gmsh.model.geo.synchronize()

    curve_tags = Int[t for (_, t) in gmsh.model.getEntities(1)]
    _setup_mesh_size!(mesh_size, dem, curve_tags)
    !isnothing(mesh_algorithm) &&
      gmsh.option.setNumber("Mesh.Algorithm", Float64(mesh_algorithm))
    recombine && gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.mesh.generate(2)
    order > 1 && gmsh.model.mesh.setOrder(order)

    # Lift every mesh node's z-coordinate by sampling the DEM.
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    n = length(node_tags)
    pts = NTuple{2,Float64}[(coords[3i-2], coords[3i-1]) for i in 1:n]
    elevations = sample_elevation(pts, dem)
    for i in 1:n
      gmsh.model.mesh.setNode(
        node_tags[i],
        [coords[3i-2], coords[3i-1], elevations[i]],
        Float64[],
      )
    end

    gmsh.write(path)

    node_tags2, _, _  = gmsh.model.mesh.getNodes()
    _, elem_tags, _   = gmsh.model.mesh.getElements(2)
    return (; nodes    = length(node_tags2),
              elements = sum(length(t) for t in elem_tags; init = 0))
  finally
    gmsh.finalize()
  end
end
