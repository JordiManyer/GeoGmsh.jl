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
(non-zero z). Interior elevation is not encoded in the `.geo` file — use
[`generate_mesh`](@ref) with a `DEMRaster` to lift interior mesh nodes.

Keyword arguments are identical to the 2D version.
"""
function write_geo(
  geoms            :: Vector{Geometry3D},
  name             :: AbstractString;
  mesh_size        :: Real               = 1.0,
  mesh_algorithm   :: Union{Int,Nothing} = nothing,
  split_components :: Bool               = false,
  verbose          :: Bool               = false,
)
  if split_components
    mkpath(name)
    n  = length(geoms)
    nd = ndigits(n)
    for (i, g) in enumerate(geoms)
      bname = (isempty(g.base.name) ? lpad(i, nd, '0') : g.base.name) * ".geo"
      _write_geo_single_3d([g], joinpath(name, bname); mesh_size, mesh_algorithm)
    end
  else
    _write_geo_single_3d(geoms, name * ".geo"; mesh_size, mesh_algorithm)
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
)
  lc = Float64(mesh_size)
  open(path, "w") do io
    _write_header(io, length(geoms), mesh_algorithm)
    pt_id = line_id = loop_id = surf_id = 1
    for g in geoms
      pt_id, line_id, loop_id, surf_id =
        _write_geometry_3d(io, g, pt_id, line_id, loop_id, surf_id, lc)
    end
  end
end

function _write_geometry_3d(io, g::Geometry3D, pt_id, line_id, loop_id, surf_id, lc)
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
  surf_id += 1
  println(io)

  return pt_id, line_id, loop_id, surf_id
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
                         mesh_size=1.0, mesh_algorithm=nothing,
                         split_components=false, verbose=false)

Generate a volumetric tetrahedral mesh by extruding the terrain surface
downward to a flat bottom plane.  Exactly one of `depth` or `z_bottom` must
be provided.

The algorithm:
1. Generate a flat 2D triangle mesh from `geoms`.
2. Sample `dem` at every mesh node to get `z_top`.
3. Create a mirror node set at the flat bottom elevation.
4. Split each triangular prism into 3 tetrahedra.
5. Tag **Top**, **Bottom**, **Sides**, and **Volume** as Gmsh physical groups
   so boundary conditions can be applied in FEM solvers.
6. Write the result as a `.msh` file.

# Arguments
- `geoms` — `Vector{Geometry3D}` (only the 2D base geometry is used).
- `dem`   — [`DEMRaster`](@ref) for elevation sampling.
- `name`  — output path **without** the `.msh` extension.

# Keyword arguments
- `depth`            — thickness below `min(z_terrain)` (same CRS units).
- `z_bottom`         — absolute bottom elevation; overrides `depth`.
- `mesh_size`        — characteristic element length (default `1.0`).
- `mesh_algorithm`   — Gmsh 2D algorithm tag (default: Gmsh's built-in default).
- `split_components` — write one `.msh` per geometry component (default `false`).
- `verbose`          — print progress (default `false`).
"""
function generate_mesh_volume(
  geoms            :: Vector{Geometry3D},
  dem              :: DEMRaster,
  name             :: AbstractString;
  depth            :: Union{Real,Nothing}              = nothing,
  z_bottom         :: Union{Real,Nothing}              = nothing,
  mesh_size        :: Union{Real,AdaptivityAlgorithm}  = 1.0,
  mesh_algorithm   :: Union{Int,Nothing}               = nothing,
  split_components :: Bool                             = false,
  verbose          :: Bool                             = false,
)
  if isnothing(depth) && isnothing(z_bottom)
    error("Either `depth` or `z_bottom` must be specified.")
  end
  depth_f   = isnothing(depth)    ? nothing : Float64(depth)
  z_bot_abs = isnothing(z_bottom) ? nothing : Float64(z_bottom)

  if split_components
    mkpath(name)
    n  = length(geoms)
    nd = ndigits(n)
    total_nodes = total_tets = 0
    for (i, g) in enumerate(geoms)
      bname = (isempty(g.base.name) ? lpad(i, nd, '0') : g.base.name) * ".msh"
      fpath = joinpath(name, bname)
      verbose && @printf("  [%*d / %d]  %-*s  ", nd, i, n, nd + 4, bname)
      stats = _generate_mesh_volume_single([g], dem, fpath;
        mesh_size, mesh_algorithm, depth = depth_f, z_bottom = z_bot_abs)
      total_nodes += stats.nodes
      total_tets  += stats.elements
      verbose && @printf("%s nodes  %s tets\n", _fmt(stats.nodes), _fmt(stats.elements))
    end
    if verbose
      println("  Total      : $(_fmt(total_nodes)) nodes  $(_fmt(total_tets)) tets")
      println("  Written    : ", name, "/")
    end
  else
    stats = _generate_mesh_volume_single(geoms, dem, name * ".msh";
      mesh_size, mesh_algorithm, depth = depth_f, z_bottom = z_bot_abs)
    if verbose
      println("  Nodes      : $(_fmt(stats.nodes))   Tets : $(_fmt(stats.elements))")
      println("  Written    : ", name, ".msh")
    end
  end
  return name
end

function _generate_mesh_volume_single(
  geoms    :: Vector{Geometry3D},
  dem      :: DEMRaster,
  path     :: AbstractString;
  mesh_size,
  mesh_algorithm,
  depth    :: Union{Float64,Nothing},
  z_bottom :: Union{Float64,Nothing},
) :: NamedTuple
  lc = _point_lc(mesh_size)

  # ── Step 1: generate flat 2D triangle mesh, capture nodes + connectivity ──
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
    tri_idx = findfirst(==(2), etypes)   # element type 2 = linear triangle
    !isnothing(tri_idx) && append!(tri_nodes_flat, enodes[tri_idx])
  finally
    gmsh.finalize()
  end

  n_nodes = length(node_tags_2d)
  n_tris  = length(tri_nodes_flat) ÷ 3
  tag_to_idx = Dict{Int,Int}(tag => i for (i, tag) in enumerate(node_tags_2d))

  # ── Step 2: sample DEM; resolve bottom elevation ───────────────────────────
  pts_2d = NTuple{2,Float64}[(coords_2d[3i-2], coords_2d[3i-1]) for i in 1:n_nodes]
  z_top  = sample_elevation(pts_2d, dem)
  z_bot  = isnothing(z_bottom) ? minimum(z_top) - depth : z_bottom

  # ── Step 3: find directed boundary edges for side-wall triangles ───────────
  # A directed edge (a → b) is a boundary edge when its reverse (b → a)
  # does not appear in any triangle (i.e. it is on the domain boundary).
  directed_edges = Set{NTuple{2,Int}}()
  for k in 1:n_tris
    t0 = tag_to_idx[tri_nodes_flat[3k-2]]
    t1 = tag_to_idx[tri_nodes_flat[3k-1]]
    t2 = tag_to_idx[tri_nodes_flat[3k  ]]
    push!(directed_edges, (t0, t1), (t1, t2), (t2, t0))
  end
  boundary_edges = NTuple{2,Int}[
    (a, b) for (a, b) in directed_edges if (b, a) ∉ directed_edges
  ]
  n_side_tris = 2 * length(boundary_edges)

  # ── Step 4: build the volumetric mesh with physical groups ─────────────────
  # Node layout: top = 1..n_nodes, bottom = n_nodes+1..2*n_nodes.
  # Element tag layout (globally unique):
  #   tets:        1          .. 3*n_tris
  #   top tris:    3*n_tris+1 .. 4*n_tris
  #   bot tris:    4*n_tris+1 .. 5*n_tris
  #   side tris:   5*n_tris+1 .. 5*n_tris+n_side_tris

  gmsh.initialize()
  try
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("GeoGmsh3DVolume")

    vol_tag  = gmsh.model.addDiscreteEntity(3)
    top_tag  = gmsh.model.addDiscreteEntity(2)
    bot_tag  = gmsh.model.addDiscreteEntity(2)
    side_tag = gmsh.model.addDiscreteEntity(2)

    # All nodes go to the volume entity.
    all_node_tags = collect(1:2*n_nodes)
    all_coords    = Vector{Float64}(undef, 6 * n_nodes)
    for i in 1:n_nodes
      x = coords_2d[3i-2];  y = coords_2d[3i-1]
      j = n_nodes + i
      all_coords[3i-2] = x;  all_coords[3i-1] = y;  all_coords[3i  ] = z_top[i]
      all_coords[3j-2] = x;  all_coords[3j-1] = y;  all_coords[3j  ] = z_bot
    end
    gmsh.model.mesh.addNodes(3, vol_tag, all_node_tags, all_coords)

    # ── Tetrahedra (prism split: 3 tets per triangle) ──────────────────────
    # For top CCW triangle (t0,t1,t2) and bottom (b0,b1,b2):
    #   Tet 1: (b0, b1, b2, t0)
    #   Tet 2: (t0, b1, b2, t2)
    #   Tet 3: (t0, t1, b1, t2)
    n_tets    = 3 * n_tris
    tet_tags  = collect(1:n_tets)
    tet_nodes = Vector{Int}(undef, 4 * n_tets)

    for k in 1:n_tris
      t0 = tag_to_idx[tri_nodes_flat[3k-2]]
      t1 = tag_to_idx[tri_nodes_flat[3k-1]]
      t2 = tag_to_idx[tri_nodes_flat[3k  ]]
      b0 = n_nodes + t0;  b1 = n_nodes + t1;  b2 = n_nodes + t2
      base = 12 * (k - 1)
      tet_nodes[base+ 1] = b0;  tet_nodes[base+ 2] = b1
      tet_nodes[base+ 3] = b2;  tet_nodes[base+ 4] = t0
      tet_nodes[base+ 5] = t0;  tet_nodes[base+ 6] = b1
      tet_nodes[base+ 7] = b2;  tet_nodes[base+ 8] = t2
      tet_nodes[base+ 9] = t0;  tet_nodes[base+10] = t1
      tet_nodes[base+11] = b1;  tet_nodes[base+12] = t2
    end
    gmsh.model.mesh.addElements(3, vol_tag, [4], [tet_tags], [tet_nodes])

    # ── Top surface triangles (same winding as 2D mesh → normal points up) ──
    top_tags  = collect(3*n_tris+1 : 4*n_tris)
    top_nodes = Vector{Int}(undef, 3 * n_tris)
    for k in 1:n_tris
      top_nodes[3k-2] = tag_to_idx[tri_nodes_flat[3k-2]]
      top_nodes[3k-1] = tag_to_idx[tri_nodes_flat[3k-1]]
      top_nodes[3k  ] = tag_to_idx[tri_nodes_flat[3k  ]]
    end
    gmsh.model.mesh.addElements(2, top_tag, [2], [top_tags], [top_nodes])

    # ── Bottom surface triangles (reversed winding → normal points down) ────
    bot_tags  = collect(4*n_tris+1 : 5*n_tris)
    bot_nodes = Vector{Int}(undef, 3 * n_tris)
    for k in 1:n_tris
      t0 = tag_to_idx[tri_nodes_flat[3k-2]]
      t1 = tag_to_idx[tri_nodes_flat[3k-1]]
      t2 = tag_to_idx[tri_nodes_flat[3k  ]]
      bot_nodes[3k-2] = n_nodes + t0
      bot_nodes[3k-1] = n_nodes + t2   # reversed: t0,t2,t1
      bot_nodes[3k  ] = n_nodes + t1
    end
    gmsh.model.mesh.addElements(2, bot_tag, [2], [bot_tags], [bot_nodes])

    # ── Side-wall triangles (two per boundary edge, outward normal) ─────────
    # For directed boundary edge (a → b, interior to the left):
    #   Tri 1: (top_a, bot_a, bot_b)
    #   Tri 2: (top_a, bot_b, top_b)
    side_tags_v  = collect(5*n_tris+1 : 5*n_tris+n_side_tris)
    side_nodes_v = Vector{Int}(undef, 3 * n_side_tris)
    for (j, (a, b)) in enumerate(boundary_edges)
      ba = n_nodes + a;  bb = n_nodes + b
      base = 6 * (j - 1)
      side_nodes_v[base+1] = a;   side_nodes_v[base+2] = ba;  side_nodes_v[base+3] = bb
      side_nodes_v[base+4] = a;   side_nodes_v[base+5] = bb;  side_nodes_v[base+6] = b
    end
    gmsh.model.mesh.addElements(2, side_tag, [2], [side_tags_v], [side_nodes_v])

    # ── Physical groups ──────────────────────────────────────────────────────
    gmsh.model.addPhysicalGroup(3, [vol_tag],  1, "Volume")
    gmsh.model.addPhysicalGroup(2, [top_tag],  2, "Top")
    gmsh.model.addPhysicalGroup(2, [bot_tag],  3, "Bottom")
    gmsh.model.addPhysicalGroup(2, [side_tag], 4, "Sides")

    gmsh.write(path)

    return (; nodes = 2 * n_nodes, elements = n_tets)
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
