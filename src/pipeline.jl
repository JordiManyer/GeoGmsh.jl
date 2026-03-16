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
- `select`            — restrict which Shapefile records are loaded.  May be
                        an `AbstractVector{Int}` of 1-based row indices, a
                        predicate `row -> Bool` on DBF attributes, or `nothing`
                        (default, load all records).  See [`list_components`](@ref)
                        to inspect available attributes.
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
- `name_fn`           — optional callable `row -> String` (same ring-level `row`
                        as `select`) that sets the filename stem for each
                        component when `split_components = true`.  When
                        `nothing` (default), files are numbered sequentially.
- `verbose`           — print progress and geometry statistics to stdout
                        (default: `true`).

Returns `output_name`.
"""
function shapefile_to_geo(
  input_path  :: AbstractString,
  output_name :: AbstractString;
  proj_method       :: Union{String,Proj.Transformation,Nothing} = "EPSG:3857",
  select                                                         = nothing,
  name_fn                                                        = nothing,
  edge_length_range :: Union{Tuple{Real,Real},Nothing}           = nothing,
  coarsen_strategy  :: Symbol                                    = :iterative,
  bbox_size         :: Union{Real,Nothing}                       = nothing,
  mesh_size         :: Real                                      = 1.0,
  mesh_algorithm    :: Union{Int,Nothing}                        = nothing,
  split_components  :: Bool                                      = false,
  verbose           :: Bool                                      = true,
)
  # --- Read ------------------------------------------------------------------
  verbose && println("Reading: ", input_path)
  geoms, source_crs = read_shapefile(input_path; select, name_fn)
  if verbose
    _print_summary(_geom_summary(geoms))
  end

  # --- Reproject -------------------------------------------------------------
  if !isnothing(proj_method)
    if verbose
      println("\nReprojecting: ", _crs_label(source_crs), " → ", _crs_label(proj_method))
    end
    geoms = project_to_meters(geoms, source_crs; target = proj_method)
    verbose && _print_bbox(_geom_summary(geoms); units = "m")
  end

  # --- Edge operations -------------------------------------------------------
  if !isnothing(edge_length_range)
    min_len, max_len = Float64.(edge_length_range)
    if verbose
      max_str = isinf(max_len) ? "∞" : @sprintf("%.3g", max_len)
      println("\nCoarsening/refining: edge ∈ [", @sprintf("%.3g", min_len), ", ", max_str, "]")
    end
    pts_before = verbose ? sum(npoints(g.exterior) for g in geoms) : 0
    n_before   = length(geoms)

    geoms = coarsen_edges(geoms, min_len; strategy = coarsen_strategy)
    geoms = refine_edges(geoms, max_len)
    n_coarsened = length(geoms)
    geoms = filter_components(geoms)

    if verbose
      pts_after = sum(npoints(g.exterior) for g in geoms)
      pct = round((1 - pts_after / pts_before) * 100; digits = 1)
      println("  Points     : $(_fmt(pts_before)) → $(_fmt(pts_after))  (−$pct%)")
      n_dropped = n_coarsened - length(geoms)
      if n_dropped > 0
        println("  Filtered   : $n_dropped degenerate component(s) dropped  →  $(length(geoms)) remaining")
      end
    end
  end

  # --- Rescale ---------------------------------------------------------------
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

  # --- Write -----------------------------------------------------------------
  verbose && println("\nWriting: ", output_name, split_components ? "/" : ".geo")
  write_geo(geoms, output_name; mesh_size, mesh_algorithm, split_components, verbose)
  return output_name
end

"""
    shapefile_to_msh(input_path, output_name; kwargs...) -> String

End-to-end pipeline: read a Shapefile, optionally reproject and resample edges,
then generate a 2-D Gmsh mesh and write a `.msh` file.  `output_name` should
be given **without** the `.msh` extension.

Accepts the same keyword arguments as [`shapefile_to_geo`](@ref)
(`proj_method`, `select`, `edge_length_range`, `coarsen_strategy`,
`bbox_size`, `mesh_size`, `mesh_algorithm`, `split_components`, `verbose`),
plus:

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
  name_fn                                                        = nothing,
  edge_length_range :: Union{Tuple{Real,Real},Nothing}           = nothing,
  coarsen_strategy  :: Symbol                                    = :iterative,
  bbox_size         :: Union{Real,Nothing}                       = nothing,
  mesh_size         :: Real                                      = 1.0,
  mesh_algorithm    :: Union{Int,Nothing}                        = nothing,
  order             :: Int                                       = 1,
  recombine         :: Bool                                      = false,
  split_components  :: Bool                                      = false,
  verbose           :: Bool                                      = true,
)
  # --- Read ------------------------------------------------------------------
  verbose && println("Reading: ", input_path)
  geoms, source_crs = read_shapefile(input_path; select, name_fn)
  if verbose
    _print_summary(_geom_summary(geoms))
  end

  # --- Reproject -------------------------------------------------------------
  if !isnothing(proj_method)
    if verbose
      println("\nReprojecting: ", _crs_label(source_crs), " → ", _crs_label(proj_method))
    end
    geoms = project_to_meters(geoms, source_crs; target = proj_method)
    verbose && _print_bbox(_geom_summary(geoms); units = "m")
  end

  # --- Edge operations -------------------------------------------------------
  if !isnothing(edge_length_range)
    min_len, max_len = Float64.(edge_length_range)
    if verbose
      max_str = isinf(max_len) ? "∞" : @sprintf("%.3g", max_len)
      println("\nCoarsening/refining: edge ∈ [", @sprintf("%.3g", min_len), ", ", max_str, "]")
    end
    pts_before = verbose ? sum(npoints(g.exterior) for g in geoms) : 0
    n_before   = length(geoms)

    geoms = coarsen_edges(geoms, min_len; strategy = coarsen_strategy)
    geoms = refine_edges(geoms, max_len)
    n_coarsened = length(geoms)
    geoms = filter_components(geoms)

    if verbose
      pts_after = sum(npoints(g.exterior) for g in geoms)
      pct = round((1 - pts_after / pts_before) * 100; digits = 1)
      println("  Points     : $(_fmt(pts_before)) → $(_fmt(pts_after))  (−$pct%)")
      n_dropped = n_coarsened - length(geoms)
      if n_dropped > 0
        println("  Filtered   : $n_dropped degenerate component(s) dropped  →  $(length(geoms)) remaining")
      end
    end
  end

  # --- Rescale ---------------------------------------------------------------
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

  # --- Mesh ------------------------------------------------------------------
  verbose && println("\nMeshing: ", length(geoms), " component(s) → ",
                     output_name, split_components ? "/" : ".msh")
  generate_mesh(geoms, output_name;
    mesh_size, mesh_algorithm, order, recombine, split_components, verbose)
  return output_name
end
