"""
Mesh adaptivity strategies for spatially-varying element sizes, and geometry
filtering utilities.

## Adaptivity extension interface

To implement a custom adaptivity strategy, subtype `AdaptivityAlgorithm` and
add methods for:

- `nominal_size(alg)  -> Float64`    — the coarsest (maximum) size the strategy
                                        may produce; used for Gmsh bounds and
                                        as the characteristic length of mesh points.
- `apply_adaptivity!(alg, dem, curve_tags) -> Union{Int,Nothing}` — build
                                        Gmsh background mesh-size fields inside an
                                        **active** Gmsh session (after geometry
                                        synchronisation).  Returns the tag of the
                                        field to set as background mesh, or `nothing`
                                        if no field is needed.

The `dem` argument is a `DEMRaster` (or `nothing` when no DEM is available, e.g.
in a purely 2-D pipeline).  `curve_tags` is an `Int` vector of all curve entity
tags in the current Gmsh model.

Both `nominal_size` and `apply_adaptivity!` are exported so that downstream
packages can extend this interface without reaching into package internals.
"""

# ============================================================================
# filter_components  (geometry filtering)
# ============================================================================

"""
    filter_components(geoms; min_points = 4) -> Vector{Geometry2D}

Remove geometrically degenerate components:

- Drop any `Geometry2D` whose exterior ring has fewer than `min_points`
  vertices.
- Strip any hole with fewer than `min_points` vertices (the exterior is kept).

The default `min_points = 4` removes 3-point (triangular) rings that survive
simplification and can cause Gmsh meshing failures.
"""
function filter_components(
  geoms      :: Vector{Geometry2D};
  min_points :: Int = 4,
) :: Vector{Geometry2D}
  out = Geometry2D[]
  for g in geoms
    npoints(g.exterior) < min_points && continue
    holes = filter(h -> npoints(h) >= min_points, g.holes)
    push!(out, Geometry2D(g.exterior, holes, g.name))
  end
  return out
end

# ============================================================================
# AdaptivityAlgorithm — abstract type and interface stubs
# ============================================================================

"""
    AdaptivityAlgorithm

Abstract base type for spatially-varying mesh-size strategies.

Pass an instance as the `mesh_size` argument of [`generate_mesh`](@ref),
[`geoms_to_msh`](@ref), [`geoms_to_msh_3d`](@ref), etc.  Passing a plain
`Real` uses a uniform characteristic length (the original behaviour).

## Implementing a custom strategy

```julia
struct MyAdaptivity <: AdaptivityAlgorithm
    coarse :: Float64
    fine   :: Float64
end

GeoGmsh.nominal_size(alg::MyAdaptivity) = alg.coarse

function GeoGmsh.apply_adaptivity!(alg::MyAdaptivity, dem, curve_tags)
    # build Gmsh fields; return the background field tag (or nothing)
    f = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(f, "F", "\$(alg.fine)")
    return f
end
```
"""
abstract type AdaptivityAlgorithm end

"""
    nominal_size(alg::AdaptivityAlgorithm) -> Float64

Return the largest (coarsest) mesh size that `alg` may produce.

This value is used as:
- the characteristic length (`lc`) passed to Gmsh point entities, and
- the `Mesh.CharacteristicLengthMax` global bound.

**Must be implemented** for every concrete `AdaptivityAlgorithm` subtype.
"""
function nominal_size end

"""
    apply_adaptivity!(alg::AdaptivityAlgorithm, dem, curve_tags) -> Union{Int,Nothing}

Build Gmsh background mesh-size fields for `alg` inside an **active** Gmsh
session (i.e. after `gmsh.initialize()` and geometry synchronisation).

Returns the Gmsh field tag to register as the background mesh, or `nothing`
if the strategy needs no field (e.g. because it degenerates to a uniform size).

# Arguments
- `dem`        — [`DEMRaster`](@ref) or `nothing` if no DEM is available
                 (2D pipeline).  Strategies that require a DEM should call
                 `error(...)` when this is `nothing`.
- `curve_tags` — `Vector{Int}` of all curve entity tags in the current model.
                 Strategies that need boundary distance should use these.

**Must be implemented** for every concrete `AdaptivityAlgorithm` subtype.
"""
function apply_adaptivity! end

# ============================================================================
# SlopeAdaptivity
# ============================================================================

"""
    SlopeAdaptivity(; mesh_size_min, mesh_size_max, percentile = 0.95)

Slope-based spatially-varying mesh size: steep terrain gets small elements,
flat terrain gets large elements.

The DEM gradient magnitude `|∇h|` is computed via finite differences on the
raster grid.  Values are normalised by the `percentile`-th percentile of the
gradient distribution (default: 95th), so outlier cliffs do not compress the
whole scale.  The resulting mesh size at each DEM pixel is:

```
t    = clamp(|∇h| / slope_at_percentile, 0, 1)
size = mesh_size_min + (mesh_size_max - mesh_size_min) * (1 − t)
```

This produces `mesh_size_max` on perfectly flat terrain and `mesh_size_min` at
slopes at or above the `percentile` threshold.

**Requires a DEM** — only usable with [`geoms_to_msh_3d`](@ref),
[`geoms_to_msh_3d_volume`](@ref), or by passing a `DEMRaster` directly to
[`generate_mesh`](@ref).

# Example
```julia
geoms_to_msh_3d("region.geojson", "dem.tif", "output";
    mesh_size  = SlopeAdaptivity(mesh_size_min = 100.0, mesh_size_max = 500.0),
    target_crs = "EPSG:32632",
)
```
"""
struct SlopeAdaptivity <: AdaptivityAlgorithm
  mesh_size_min :: Float64
  mesh_size_max :: Float64
  percentile    :: Float64
end

SlopeAdaptivity(;
  mesh_size_min :: Real,
  mesh_size_max :: Real,
  percentile    :: Real = 0.95,
) = SlopeAdaptivity(Float64(mesh_size_min), Float64(mesh_size_max), Float64(percentile))

nominal_size(alg::SlopeAdaptivity) = alg.mesh_size_max

function apply_adaptivity!(alg::SlopeAdaptivity, ::Nothing, _)
  error("SlopeAdaptivity requires a DEM. " *
        "Use geoms_to_msh_3d / geoms_to_msh_3d_volume, or pass a DEMRaster " *
        "directly to generate_mesh.")
end

function apply_adaptivity!(alg::SlopeAdaptivity, dem::DEMRaster, _curve_tags)
  gt    = dem.transform
  data  = dem.data
  ncols, nrows = size(data)
  dx = abs(gt[2])
  dy = abs(gt[6])

  # Gradient magnitude via central differences (forward/backward at borders).
  grad = Matrix{Float64}(undef, ncols, nrows)
  for r in 1:nrows, c in 1:ncols
    gx = c == 1     ? (data[c+1, r] - data[c,   r]) / dx :
         c == ncols ? (data[c,   r] - data[c-1, r]) / dx :
                      (data[c+1, r] - data[c-1, r]) / (2dx)
    gy = r == 1     ? (data[c, r+1] - data[c, r  ]) / dy :
         r == nrows ? (data[c, r  ] - data[c, r-1]) / dy :
                      (data[c, r+1] - data[c, r-1]) / (2dy)
    grad[c, r] = sqrt(gx^2 + gy^2)
  end

  # Normalise by the requested percentile of the gradient distribution.
  sorted  = sort!(vec(copy(grad)))
  idx     = clamp(round(Int, alg.percentile * length(sorted)), 1, length(sorted))
  g_scale = max(sorted[idx], 1e-10)

  smin, smax = alg.mesh_size_min, alg.mesh_size_max

  # Build a Gmsh PostView: one scalar value per DEM pixel (SP = scalar point).
  n_pts   = ncols * nrows
  sp_data = Vector{Float64}(undef, 4 * n_pts)
  k = 0
  for r in 1:nrows, c in 1:ncols
    x = gt[1] + (c - 0.5) * gt[2]   # pixel-centre world x
    y = gt[4] + (r - 0.5) * gt[6]   # pixel-centre world y (gt[6] < 0 ⟹ north-up)
    t = clamp(grad[c, r] / g_scale, 0.0, 1.0)
    s = smin + (smax - smin) * (1.0 - t)
    sp_data[k+1] = x
    sp_data[k+2] = y
    sp_data[k+3] = 0.0
    sp_data[k+4] = s
    k += 4
  end

  view_tag = gmsh.view.add("GeoGmsh_slope")
  gmsh.view.addListData(view_tag, "SP", n_pts, sp_data)

  f = gmsh.model.mesh.field.add("PostView")
  gmsh.model.mesh.field.setNumber(f, "ViewTag", Float64(view_tag))
  return f
end

# ============================================================================
# BoundaryLayerAdaptivity
# ============================================================================

"""
    BoundaryLayerAdaptivity(; mesh_size, mesh_size_min, width)

Boundary-layer mesh refinement: small elements right at polygon edges,
smoothly coarsening to `mesh_size` over a distance of `width` (in the mesh
CRS units, typically metres for a UTM projection).

Uses Gmsh's `Distance` + `Threshold` fields:

```
size(d) = mesh_size_min   for d ≤ 0
         ↗ (linear ramp)
        mesh_size          for d ≥ width
```

Works in both 2D and 3D pipelines (curve tags are available in both).

# Example
```julia
geoms_to_msh("region.geojson", "output";
    mesh_size = BoundaryLayerAdaptivity(
        mesh_size     = 500.0,
        mesh_size_min =  50.0,
        width         = 2_000.0,
    ),
)
```
"""
struct BoundaryLayerAdaptivity <: AdaptivityAlgorithm
  mesh_size     :: Float64
  mesh_size_min :: Float64
  width         :: Float64
end

BoundaryLayerAdaptivity(;
  mesh_size     :: Real,
  mesh_size_min :: Real,
  width         :: Real,
) = BoundaryLayerAdaptivity(Float64(mesh_size), Float64(mesh_size_min), Float64(width))

nominal_size(alg::BoundaryLayerAdaptivity) = alg.mesh_size

function apply_adaptivity!(alg::BoundaryLayerAdaptivity, _dem, curve_tags)
  isempty(curve_tags) && return nothing

  d = gmsh.model.mesh.field.add("Distance")
  gmsh.model.mesh.field.setNumbers(d, "CurvesList", Float64.(curve_tags))

  t = gmsh.model.mesh.field.add("Threshold")
  gmsh.model.mesh.field.setNumber(t, "InField",  Float64(d))
  gmsh.model.mesh.field.setNumber(t, "SizeMin",  alg.mesh_size_min)
  gmsh.model.mesh.field.setNumber(t, "SizeMax",  alg.mesh_size)
  gmsh.model.mesh.field.setNumber(t, "DistMin",  0.0)
  gmsh.model.mesh.field.setNumber(t, "DistMax",  alg.width)
  return t
end

# ============================================================================
# ComposedAdaptivity
# ============================================================================

"""
    ComposedAdaptivity(alg1, alg2, ...)
    ComposedAdaptivity(algorithms::Vector{AdaptivityAlgorithm})

Combine multiple `AdaptivityAlgorithm`s using element-wise minimum (the finest
size wins at every point in the domain).  Implemented with Gmsh's `Min` field.

# Example — slope refinement + boundary layer together
```julia
geoms_to_msh_3d("region.geojson", "dem.tif", "output";
    mesh_size = ComposedAdaptivity(
        SlopeAdaptivity(mesh_size_min = 100.0, mesh_size_max = 500.0),
        BoundaryLayerAdaptivity(mesh_size = 500.0, mesh_size_min = 50.0, width = 1_000.0),
    ),
    target_crs = "EPSG:32632",
)
```
"""
struct ComposedAdaptivity <: AdaptivityAlgorithm
  algorithms :: Vector{AdaptivityAlgorithm}
end

ComposedAdaptivity(algs::AdaptivityAlgorithm...) = ComposedAdaptivity(collect(algs))

nominal_size(alg::ComposedAdaptivity) = maximum(nominal_size, alg.algorithms)

function apply_adaptivity!(alg::ComposedAdaptivity, dem, curve_tags)
  sub_tags = Int[]
  for a in alg.algorithms
    tag = apply_adaptivity!(a, dem, curve_tags)
    isnothing(tag) || push!(sub_tags, tag)
  end
  isempty(sub_tags) && return nothing
  length(sub_tags) == 1 && return only(sub_tags)

  f = gmsh.model.mesh.field.add("Min")
  gmsh.model.mesh.field.setNumbers(f, "FieldsList", Float64.(sub_tags))
  return f
end
