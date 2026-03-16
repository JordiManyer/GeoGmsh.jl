module ShapefileToGmsh

using Shapefile
using Printf

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
- `proj_method`       — how to reproject coordinates.  May be a
                        `ProjectionMethod` struct (`Equirectangular()`,
                        `Mercator()`, `WebMercator()`), a convenience symbol
                        (`:equirectangular`, `:mercator`, `:web_mercator`), or
                        `nothing` to skip reprojection entirely (use when the
                        shapefile is already in metres).  Default:
                        `:equirectangular`.
- `edge_length_range` — `(min, max)` edge length bounds in the same units as
                        the (possibly reprojected) coordinates.  Edges shorter
                        than `min` are coarsened; edges longer than `max` are
                        refined.  `nothing` skips both steps.  Default:
                        `nothing`.
- `coarsen_strategy`  — `:single` or `:iterative` (default: `:iterative`).
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
  proj_method      :: Union{ProjectionMethod,Symbol,Nothing} = :equirectangular,
  edge_length_range :: Union{Tuple{Real,Real},Nothing}       = nothing,
  coarsen_strategy :: Symbol                                 = :iterative,
  mesh_size        :: Real                                   = 1.0,
  mesh_algorithm   :: Union{Int,Nothing}                     = nothing,
  split_components :: Bool                                   = false,
)
  geoms, crs = read_shapefile(input_path)

  if !isnothing(proj_method)
    geoms = project_to_meters(geoms, crs; method=proj_method)
  end

  if !isnothing(edge_length_range)
    min_len, max_len = Float64.(edge_length_range)
    geoms = coarsen_edges(geoms, min_len; strategy=coarsen_strategy)
    geoms = refine_edges(geoms, max_len)
  end

  write_geo(geoms, output_name;
    mesh_size, mesh_algorithm, split_components)
  return output_name
end

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

export ShapeGeometry, Contour, npoints, nedges

export read_shapefile

export ProjectionMethod, Equirectangular, Mercator, WebMercator
export project_to_meters

export coarsen_edges, refine_edges

export write_geo

export shapefile_to_geo

end
