"""
    geoms_to_geo(geom, output_name; kwargs...) -> String

Convert any GeoInterface-compatible geometry (or a `DataFrame` with a geometry
column) to a Gmsh `.geo` script.  `output_name` should be given **without**
the `.geo` extension.

The pipeline runs these steps before writing:
1. Reproject (`GeometryOps.reproject`) — if `target_crs` is set
2. Simplify (`GeometryOps.simplify` with [`MinEdgeLength`](@ref)) — if `simplify_tol` is set
3. Segmentize (`GeometryOps.segmentize`) — if `max_edge_length` is set
4. Ingest ([`ingest`](@ref)) — convert to Gmsh-ready internal representation
5. Filter ([`filter_components`](@ref)) — drop degenerate rings
6. Rescale ([`rescale`](@ref)) — if `bbox_size` is set

# Keyword arguments
- `target_crs`       — destination CRS string (e.g. `"EPSG:3857"`) or `nothing`
                       to skip reprojection.  Default: `"EPSG:3857"`.
- `select`           — when `geom` is a file path, a predicate `row -> Bool`
                       passed to [`read_geodata`](@ref).  Ignored otherwise.
- `simplify_tol`     — minimum edge length after simplification (in the
                       coordinate units after reprojection).  `nothing` skips.
- `max_edge_length`  — maximum edge length for segmentization.  `nothing` skips.
- `bbox_size`        — rescale so the largest bounding-box dimension equals
                       this value, origin at (0, 0).  `nothing` to skip.
- `mesh_size`        — characteristic element length (default: `1.0`).
- `mesh_algorithm`   — Gmsh algorithm tag (default: `nothing`).
- `split_components` — write one `.geo` file per component (default: `false`).
- `verbose`          — print progress (default: `true`).

Returns `output_name`.
"""
function geoms_to_geo(
  geom,
  output_name :: AbstractString;
  target_crs        :: Union{String,Nothing} = "EPSG:3857",
  select                                     = nothing,
  simplify_tol      :: Union{Real,Nothing}   = nothing,
  max_edge_length   :: Union{Real,Nothing}   = nothing,
  bbox_size         :: Union{Real,Nothing}   = nothing,
  mesh_size         :: Real                  = 1.0,
  mesh_algorithm    :: Union{Int,Nothing}    = nothing,
  split_components  :: Bool                  = false,
  verbose           :: Bool                  = true,
)
  geom, source_crs = _load(geom; select, verbose)
  geoms = _run_pipeline(geom, source_crs;
    target_crs, simplify_tol, max_edge_length, bbox_size, verbose)
  verbose && println("\nWriting: ", output_name, split_components ? "/" : ".geo")
  write_geo(geoms, output_name; mesh_size, mesh_algorithm, split_components, verbose)
  return output_name
end

"""
    geoms_to_msh(geom, output_name; kwargs...) -> String

Same as [`geoms_to_geo`](@ref) but generates a 2-D Gmsh mesh (`.msh` file).

Additional keyword arguments (beyond those of `geoms_to_geo`):
- `order`     — element order: 1 = linear (default), 2 = quadratic.
- `recombine` — recombine triangles into quadrilaterals (default: `false`).
"""
function geoms_to_msh(
  geom,
  output_name :: AbstractString;
  target_crs        :: Union{String,Nothing} = "EPSG:3857",
  select                                     = nothing,
  simplify_tol      :: Union{Real,Nothing}   = nothing,
  max_edge_length   :: Union{Real,Nothing}   = nothing,
  bbox_size         :: Union{Real,Nothing}   = nothing,
  mesh_size         :: Real                  = 1.0,
  mesh_algorithm    :: Union{Int,Nothing}    = nothing,
  order             :: Int                   = 1,
  recombine         :: Bool                  = false,
  split_components  :: Bool                  = false,
  verbose           :: Bool                  = true,
)
  geom, source_crs = _load(geom; select, verbose)
  geoms = _run_pipeline(geom, source_crs;
    target_crs, simplify_tol, max_edge_length, bbox_size, verbose)
  verbose && println("\nMeshing: ", length(geoms), " component(s) → ",
                     output_name, split_components ? "/" : ".msh")
  generate_mesh(geoms, output_name;
    mesh_size, mesh_algorithm, order, recombine, split_components, verbose)
  return output_name
end

# ---------------------------------------------------------------------------
# 3D terrain pipeline
# ---------------------------------------------------------------------------

"""
    geoms_to_geo_3d(geom, dem_path, output_name; kwargs...) -> String

Terrain-aware variant of [`geoms_to_geo`](@ref): runs the standard 2D
pipeline, reads the DEM at `dem_path`, lifts boundary points to terrain
elevation, and writes a `.geo` file with 3D point coordinates.

The DEM must be in the same CRS as the (reprojected) vector data.

All 2D keyword arguments are supported. Additional argument:
- `nodata_fill` — elevation substituted for out-of-bounds or nodata cells
                  (default `0.0`).
"""
function geoms_to_geo_3d(
  geom,
  dem_path     :: AbstractString,
  output_name  :: AbstractString;
  nodata_fill  :: Real = 0.0,
  verbose      :: Bool = true,
  kwargs...,
)
  geom_raw, source_crs = _load(geom; select = get(kwargs, :select, nothing), verbose)
  geoms2d = _run_pipeline(geom_raw, source_crs;
    target_crs      = get(kwargs, :target_crs,      "EPSG:3857"),
    simplify_tol    = get(kwargs, :simplify_tol,    nothing),
    max_edge_length = get(kwargs, :max_edge_length,  nothing),
    bbox_size       = get(kwargs, :bbox_size,        nothing),
    verbose,
  )
  verbose && println("\nReading DEM: ", dem_path)
  dem    = read_dem(dem_path)
  geoms3d = lift_to_3d(geoms2d, dem; nodata_fill = Float64(nodata_fill))
  verbose && println("\nWriting: ", output_name, ".geo")
  write_geo(geoms3d, output_name;
    mesh_size      = get(kwargs, :mesh_size,      1.0),
    mesh_algorithm = get(kwargs, :mesh_algorithm, nothing),
    split_components = get(kwargs, :split_components, false),
    verbose,
  )
  return output_name
end

"""
    geoms_to_msh_3d(geom, dem_path, output_name; kwargs...) -> String

Terrain-aware variant of [`geoms_to_msh`](@ref): runs the standard 2D
pipeline to generate a flat mesh, then lifts every mesh node's z-coordinate
by sampling `dem_path` at its (x, y) position.

The DEM must be in the same CRS as the (reprojected) vector data.

All 2D keyword arguments are supported. Additional argument:
- `nodata_fill` — elevation substituted for out-of-bounds or nodata cells
                  (default `0.0`).
"""
function geoms_to_msh_3d(
  geom,
  dem_path     :: AbstractString,
  output_name  :: AbstractString;
  nodata_fill  :: Real = 0.0,
  verbose      :: Bool = true,
  kwargs...,
)
  geom_raw, source_crs = _load(geom; select = get(kwargs, :select, nothing), verbose)
  geoms2d = _run_pipeline(geom_raw, source_crs;
    target_crs      = get(kwargs, :target_crs,      "EPSG:3857"),
    simplify_tol    = get(kwargs, :simplify_tol,    nothing),
    max_edge_length = get(kwargs, :max_edge_length,  nothing),
    bbox_size       = get(kwargs, :bbox_size,        nothing),
    verbose,
  )
  verbose && println("\nReading DEM: ", dem_path)
  dem     = read_dem(dem_path)
  geoms3d = lift_to_3d(geoms2d, dem; nodata_fill = Float64(nodata_fill))
  verbose && println("\nMeshing (3D): ", length(geoms3d), " component(s) → ",
                     output_name, ".msh")
  generate_mesh(geoms3d, dem, output_name;
    mesh_size      = get(kwargs, :mesh_size,      1.0),
    mesh_algorithm = get(kwargs, :mesh_algorithm, nothing),
    order          = get(kwargs, :order,          1),
    recombine      = get(kwargs, :recombine,      false),
    split_components = get(kwargs, :split_components, false),
    verbose,
  )
  return output_name
end

# ---------------------------------------------------------------------------
# Backward-compatible wrappers
# ---------------------------------------------------------------------------

"""
    shapefile_to_geo(path, output_name; kwargs...) -> String

Backward-compatible wrapper: read a geospatial file and produce a `.geo` file.
All keyword arguments are forwarded to [`geoms_to_geo`](@ref).
"""
shapefile_to_geo(path::AbstractString, out::AbstractString; kwargs...) =
  geoms_to_geo(read_geodata(path), out; kwargs...)

"""
    shapefile_to_msh(path, output_name; kwargs...) -> String

Backward-compatible wrapper: read a geospatial file and produce a `.msh` file.
All keyword arguments are forwarded to [`geoms_to_msh`](@ref).
"""
shapefile_to_msh(path::AbstractString, out::AbstractString; kwargs...) =
  geoms_to_msh(read_geodata(path), out; kwargs...)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Load from any supported input type.
# Returns (raw_geom_or_df, source_crs_wkt::Union{String,Nothing}).
function _load(path::AbstractString; select, verbose)
  verbose && println("Reading: ", path)
  df = read_geodata(path; select)
  _load(df; select = nothing, verbose = false)
end

function _load(df::DataFrames.AbstractDataFrame; select, verbose)
  df2 = isnothing(select) ? df : filter(select, df)
  verbose && println("Reading: DataFrame ($(nrow(df2)) rows)")
  crs        = GI.crs(df2)
  source_crs = _crs_to_wkt(crs)
  return df2, source_crs
end

function _load(geom; select, verbose)
  verbose && println("Ingesting: $(typeof(geom))")
  return geom, nothing   # no CRS available from raw GI geometries
end

# Apply f to each geometry in a DataFrame's geometry column.
function _apply_to_geoms(df::DataFrames.AbstractDataFrame, f)
  col    = first(GI.geometrycolumns(df))
  result = copy(df)
  result[!, col] = map(g -> isnothing(g) ? g : f(g), result[!, col])
  return result
end

# Apply f directly to a raw GI geometry.
_apply_to_geoms(geom, f) = f(geom)

# Run the full pre-Gmsh pipeline on raw GI geometry / DataFrame.
function _run_pipeline(
  geom, source_crs;
  target_crs, simplify_tol, max_edge_length, bbox_size, verbose,
)
  # --- Reproject (on raw GI geometry, before ingest) ---
  if !isnothing(target_crs)
    if isnothing(source_crs)
      @warn "No CRS information found; assuming EPSG:4326 (geographic degrees)."
      source_crs = "EPSG:4326"
    end
    verbose && println("\nReprojecting: ", _crs_label(source_crs), " → ", _crs_label(target_crs))
    geom = _apply_to_geoms(geom, g -> GO.reproject(g; source_crs, target_crs))
  end

  # --- Simplify (on raw GI geometry, before ingest) ---
  if !isnothing(simplify_tol)
    alg = MinEdgeLength(tol = Float64(simplify_tol))
    verbose && println("\nSimplifying: min edge length = $simplify_tol")
    geom = _apply_to_geoms(geom, g -> GO.simplify(alg, g))
  end

  # --- Segmentize (on raw GI geometry, before ingest) ---
  if !isnothing(max_edge_length)
    md = Float64(max_edge_length)
    verbose && println("\nSegmentizing: max edge length = $md")
    geom = _apply_to_geoms(geom, g -> GO.segmentize(g; max_distance = md))
  end

  # --- Ingest → Gmsh-ready structs (last step before Gmsh) ---
  geoms = ingest(geom)
  geoms = filter_components(geoms)
  verbose && _print_summary(_geom_summary(geoms))

  # --- Rescale (post-ingest; Gmsh-specific normalisation) ---
  if !isnothing(bbox_size)
    verbose && println("\nRescaling: L = $bbox_size")
    xmin0, xmax0, ymin0, ymax0 = _global_bbox(geoms)
    scale = Float64(bbox_size) / max(xmax0 - xmin0, ymax0 - ymin0)
    geoms = rescale(geoms, bbox_size)
    if verbose
      @printf("  Scale      : %.4g\n", scale)
      _print_bbox(_geom_summary(geoms))
    end
  end

  return geoms
end
