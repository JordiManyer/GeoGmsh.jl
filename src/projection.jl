"""
Coordinate projection utilities.

All projections convert geographic coordinates (lon, lat in decimal degrees)
to planar Cartesian coordinates in metres.

## Projection method types

Three concrete subtypes of `ProjectionMethod` are provided:

| Type                    | Description                                            |
|-------------------------|--------------------------------------------------------|
| `Equirectangular`       | Simple plate carrée. Fast; accurate for small regions. |
| `Mercator`              | Conformal. Preserves angles; distorts area at high lat.|
| `WebMercator`           | Identical formula to Mercator; labelled EPSG:3857.     |

For convenience, the corresponding symbols `:equirectangular`, `:mercator`,
and `:web_mercator` are also accepted wherever a `ProjectionMethod` is expected.

## Equirectangular

    x = R · cos(φ₀) · λ
    y = R · φ

where λ and φ are longitude and latitude in radians and φ₀ is the reference
(central) latitude. Constructed as `Equirectangular()` (auto-centroid) or
`Equirectangular(lat_ref = φ₀)`.

## Mercator / Web Mercator

    x = R · λ
    y = R · ln( tan(π/4 + φ/2) )

Mercator is ill-defined at ±90°; input latitudes are clamped to ±85.051°.
"""

# ============================================================================
# Public types
# ============================================================================

abstract type ProjectionMethod end

"""
    Equirectangular(; lat_ref = nothing)

Equirectangular (plate carrée) projection.
`lat_ref` is the reference latitude in degrees; when `nothing` it is
auto-computed as the centroid latitude of the input geometries.
"""
struct Equirectangular <: ProjectionMethod
  lat_ref::Union{Float64,Nothing}
end
Equirectangular(; lat_ref = nothing) = Equirectangular(lat_ref)

"""    Mercator()

Mercator conformal projection."""
struct Mercator <: ProjectionMethod end

"""    WebMercator()

Web Mercator projection (EPSG:3857). Identical formula to `Mercator`."""
struct WebMercator <: ProjectionMethod end

# ============================================================================
# Public API
# ============================================================================

const _EARTH_RADIUS = 6_371_000.0   # metres (mean radius)

"""
    project_to_meters(geoms, crs = :degrees; method = Equirectangular())
    -> Vector{ShapeGeometry}

Reproject `geoms` from degrees to metres.

`method` may be a `ProjectionMethod` struct or one of the convenience symbols
`:equirectangular`, `:mercator`, `:web_mercator`.

If `crs == :meters` the geometries are returned unchanged (identity).
If `crs == :unknown` a warning is issued and the projection proceeds anyway.
"""
function project_to_meters(
  geoms  :: Vector{ShapeGeometry},
  crs    :: Symbol = :degrees;
  method :: Union{ProjectionMethod,Symbol} = Equirectangular(),
) :: Vector{ShapeGeometry}
  crs == :meters && return geoms
  crs == :unknown &&
    @warn "CRS is unknown; assuming geographic (degrees) and projecting anyway."
  m       = _to_proj_method(method)
  project = _build_projector(geoms, m)
  return [_project_geometry(g, project) for g in geoms]
end

# ============================================================================
# Internal helpers
# ============================================================================

# Convert a Symbol shorthand to the corresponding struct.
_to_proj_method(m::ProjectionMethod) = m
function _to_proj_method(s::Symbol)
  if s == :equirectangular; return Equirectangular()
  elseif s == :mercator;    return Mercator()
  elseif s == :web_mercator; return WebMercator()
  else
    throw(ArgumentError(
      "Unknown projection method: $s. " *
      "Choose from :equirectangular, :mercator, :web_mercator."
    ))
  end
end

function _build_projector(geoms, m::Equirectangular)
  φ₀    = isnothing(m.lat_ref) ? _centroid_lat(geoms) : Float64(m.lat_ref)
  cosφ₀ = cos(deg2rad(φ₀))
  return pt -> (
    _EARTH_RADIUS * cosφ₀ * deg2rad(pt[1]),
    _EARTH_RADIUS * deg2rad(pt[2]),
  )
end

function _build_projector(::Any, ::Union{Mercator,WebMercator})
  return pt -> _mercator_pt(pt)
end

function _mercator_pt(pt::NTuple{2,Float64})
  lon, lat = pt
  lat = clamp(lat, -85.051129, 85.051129)   # avoid ±Inf at poles
  x = _EARTH_RADIUS * deg2rad(lon)
  y = _EARTH_RADIUS * log(tan(π/4 + deg2rad(lat)/2))
  return (x, y)
end

function _project_contour(c::Contour, f)
  Contour(map(f, c.points), c.closed)
end

function _project_geometry(g::ShapeGeometry, f)
  ShapeGeometry(
    _project_contour(g.exterior, f),
    [_project_contour(h, f) for h in g.holes],
  )
end

# Centroid latitude: mean of all y-coordinates across all geometries.
function _centroid_lat(geoms::Vector{ShapeGeometry})
  total = 0.0
  count = 0
  for g in geoms
    for pt in g.exterior.points
      total += pt[2]
      count += 1
    end
  end
  count == 0 && return 0.0
  return total / count
end
