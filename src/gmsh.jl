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

Because `ShapeGeometry` already normalises ring orientation (see `shapefiles.jl`),
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

When `split_components = true`, one `.geo` file per `ShapeGeometry` is written
into a directory named `name/`.  Files inside are named `1.geo`, `2.geo`, … .

# Keyword arguments
- `mesh_size`        — characteristic element length assigned to every point.
- `mesh_algorithm`   — optional integer passed as `Mesh.Algorithm` (e.g. 5 =
                       Delaunay, 6 = Frontal-Delaunay, 8 = Frontal-Quad).
                       Omitted from the file when `nothing`.
- `split_components` — write one file per geometry component (default `false`).
"""
function write_geo(
  geoms            :: Vector{ShapeGeometry},
  name             :: AbstractString;
  mesh_size        :: Real               = 1.0,
  mesh_algorithm   :: Union{Int,Nothing} = nothing,
  split_components :: Bool               = false,
)
  if split_components
    _write_geo_split(geoms, name; mesh_size, mesh_algorithm)
  else
    _write_geo_single(geoms, name * ".geo"; mesh_size, mesh_algorithm)
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

When `split_components = true`, one `.msh` file per `ShapeGeometry` is written
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
  geoms            :: Vector{ShapeGeometry},
  name             :: AbstractString;
  mesh_size        :: Real               = 1.0,
  mesh_algorithm   :: Union{Int,Nothing} = nothing,
  order            :: Int                = 1,
  recombine        :: Bool               = false,
  split_components :: Bool               = false,
)
  if split_components
    _generate_mesh_split(geoms, name; mesh_size, mesh_algorithm, order, recombine)
  else
    _generate_mesh_single(geoms, name * ".msh"; mesh_size, mesh_algorithm, order, recombine)
  end
  return name
end

# ============================================================================
# write_geo internals
# ============================================================================

function _write_geo_single(
  geoms          :: Vector{ShapeGeometry},
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
  geoms          :: Vector{ShapeGeometry},
  name           :: AbstractString;
  mesh_size      :: Real,
  mesh_algorithm,
)
  mkpath(name)
  nd = ndigits(length(geoms))
  for (i, g) in enumerate(geoms)
    fname = joinpath(name, lpad(i, nd, '0') * ".geo")
    _write_geo_single([g], fname; mesh_size, mesh_algorithm)
  end
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
  geoms          :: Vector{ShapeGeometry},
  path           :: AbstractString;
  mesh_size,
  mesh_algorithm,
  order,
  recombine,
)
  lc = Float64(mesh_size)
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

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
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
  finally
    gmsh.finalize()
  end
end

function _generate_mesh_split(
  geoms          :: Vector{ShapeGeometry},
  name           :: AbstractString;
  mesh_size,
  mesh_algorithm,
  order,
  recombine,
)
  mkpath(name)
  nd = ndigits(length(geoms))
  for (i, g) in enumerate(geoms)
    fname = joinpath(name, lpad(i, nd, '0') * ".msh")
    _generate_mesh_single([g], fname; mesh_size, mesh_algorithm, order, recombine)
  end
end

# Add one ShapeGeometry to the active Gmsh model; return updated counters.
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
