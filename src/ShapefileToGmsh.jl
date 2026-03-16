module ShapefileToGmsh

using Shapefile
using Printf
import Proj
import Tables

include("geometry.jl")
include("shapefiles.jl")
include("projection.jl")
include("adaptivity.jl")
include("gmsh.jl")

# ---------------------------------------------------------------------------
# Convenience pipeline
# ---------------------------------------------------------------------------

"""
    shapefile_to_geo(input_path, output_name; kwargs...) -> String

End-to-end pipeline: read a Shapefile, optionally reproject and resample edges,
then write a Gmsh `.geo` file.  `output_name` should be given **without** the
`.geo` extension.

# Keyword arguments
- `proj_method`       — target CRS string (e.g. `"EPSG:3857"`), a pre-built
                        `Proj.Transformation`, or `nothing` to skip reprojection
                        (use when the shapefile is already in the desired units).
                        Default: `"EPSG:3857"`.
- `edge_length_range` — `(min, max)` edge length bounds in the (possibly
                        reprojected) coordinate units.  Edges shorter than `min`
                        are coarsened; edges longer than `max` are refined.
                        `nothing` skips both steps.  Default: `nothing`.
- `coarsen_strategy`  — `:single` or `:iterative` (default: `:iterative`).
- `bbox_size`         — if set, rescale all geometries so the largest bounding-
                        box dimension equals this value, with the origin at
                        (0, 0).  Applied after edge operations.  Default:
                        `nothing`.
- `mesh_size`         — characteristic element length in the `.geo` file
                        (default: `1.0`).
- `mesh_algorithm`    — Gmsh meshing algorithm tag written into the `.geo` file
                        (default: `nothing`, let Gmsh decide).
- `split_components`  — if `true`, write one `.geo` file per geometry component
                        into a directory named `output_name` (default: `false`).

Returns `output_name`.
"""
function shapefile_to_geo(
  input_path  :: AbstractString,
  output_name :: AbstractString;
  proj_method       :: Union{String,Proj.Transformation,Nothing} = "EPSG:3857",
  select                                                         = nothing,
  edge_length_range :: Union{Tuple{Real,Real},Nothing}           = nothing,
  coarsen_strategy  :: Symbol                                    = :iterative,
  bbox_size         :: Union{Real,Nothing}                       = nothing,
  mesh_size         :: Real                                      = 1.0,
  mesh_algorithm    :: Union{Int,Nothing}                        = nothing,
  split_components  :: Bool                                      = false,
)
  geoms, source_crs = read_shapefile(input_path; select)

  if !isnothing(proj_method)
    geoms = project_to_meters(geoms, source_crs; target = proj_method)
  end

  if !isnothing(edge_length_range)
    min_len, max_len = Float64.(edge_length_range)
    geoms = coarsen_edges(geoms, min_len; strategy = coarsen_strategy)
    geoms = refine_edges(geoms, max_len)
  end

  if !isnothing(bbox_size)
    geoms = rescale(geoms, bbox_size)
  end

  write_geo(geoms, output_name; mesh_size, mesh_algorithm, split_components)
  return output_name
end

"""
    shapefile_to_msh(input_path, output_name; kwargs...) -> String

End-to-end pipeline: read a Shapefile, optionally reproject and resample edges,
then generate a 2-D Gmsh mesh and write a `.msh` file.  `output_name` should
be given **without** the `.msh` extension.

Accepts the same geometry kwargs as `shapefile_to_geo`, plus:

# Additional keyword arguments
- `order`     — element order: 1 = linear (default), 2 = quadratic.
- `recombine` — recombine triangles into quadrilaterals (default `false`).

Returns `output_name`.
"""
function shapefile_to_msh(
  input_path  :: AbstractString,
  output_name :: AbstractString;
  proj_method       :: Union{String,Proj.Transformation,Nothing} = "EPSG:3857",
  select                                                         = nothing,
  edge_length_range :: Union{Tuple{Real,Real},Nothing}           = nothing,
  coarsen_strategy  :: Symbol                                    = :iterative,
  bbox_size         :: Union{Real,Nothing}                       = nothing,
  mesh_size         :: Real                                      = 1.0,
  mesh_algorithm    :: Union{Int,Nothing}                        = nothing,
  order             :: Int                                       = 1,
  recombine         :: Bool                                      = false,
  split_components  :: Bool                                      = false,
)
  geoms, source_crs = read_shapefile(input_path; select)

  if !isnothing(proj_method)
    geoms = project_to_meters(geoms, source_crs; target = proj_method)
  end

  if !isnothing(edge_length_range)
    min_len, max_len = Float64.(edge_length_range)
    geoms = coarsen_edges(geoms, min_len; strategy = coarsen_strategy)
    geoms = refine_edges(geoms, max_len)
  end

  if !isnothing(bbox_size)
    geoms = rescale(geoms, bbox_size)
  end

  generate_mesh(geoms, output_name;
    mesh_size, mesh_algorithm, order, recombine, split_components)
  return output_name
end

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

export ShapeGeometry, Contour, npoints, nedges

export read_shapefile, list_components

export project_to_meters, rescale

export coarsen_edges, refine_edges

export write_geo
export generate_mesh

export shapefile_to_geo
export shapefile_to_msh

end
