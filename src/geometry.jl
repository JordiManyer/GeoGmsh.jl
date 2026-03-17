"""
    Contour

A closed or open sequence of 2D points. Closed contours represent polygon
rings; open contours represent polylines.

Points are stored as `NTuple{2,Float64}` (x, y).
"""
struct Contour
  points::Vector{NTuple{2,Float64}}
  closed::Bool
end

"""
    Geometry2D

A single polygon: an exterior ring plus zero or more hole rings.
All rings are stored as `Contour` with `closed = true`.
Exterior rings are oriented CCW (positive signed area);
hole rings are oriented CW (negative signed area).
"""
struct Geometry2D
  exterior::Contour
  holes::Vector{Contour}
  name::String
end

Geometry2D(exterior::Contour, holes::Vector{Contour}) =
  Geometry2D(exterior, holes, "")

# ---------------------------------------------------------------------------
# GeoInterface traits
# ---------------------------------------------------------------------------

import GeoInterface as GI

GI.isgeometry(::Type{Geometry2D}) = true
GI.geomtrait(::Geometry2D)        = GI.PolygonTrait()
GI.nhole(::GI.PolygonTrait, g::Geometry2D)        = length(g.holes)
GI.getexterior(::GI.PolygonTrait, g::Geometry2D)  = g.exterior
GI.gethole(::GI.PolygonTrait, g::Geometry2D, i)   = g.holes[i]

GI.isgeometry(::Type{Contour}) = true
GI.geomtrait(::Contour)        = GI.LinearRingTrait()
GI.npoint(::GI.LinearRingTrait, c::Contour)      = length(c.points)
GI.getpoint(::GI.LinearRingTrait, c::Contour, i) = c.points[i]
# LinearRingTrait is closed by definition; no need to implement isclosed.
# NTuple{2,Float64} points already have default GI.PointTrait + GI.x/GI.y from GeoInterface.

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

"""Number of points in a contour."""
npoints(c::Contour) = length(c.points)

"""Number of edges in a contour."""
nedges(c::Contour) = c.closed ? npoints(c) : npoints(c) - 1

"""Total number of points across exterior + holes."""
function npoints(g::Geometry2D)
  n = npoints(g.exterior)
  for h in g.holes
    n += npoints(h)
  end
  return n
end

# ---------------------------------------------------------------------------
# Geometry3D
# ---------------------------------------------------------------------------

"""
    Geometry3D

A terrain-aware 2+1D polygon: a `Geometry2D` base (2D boundary topology)
combined with per-point elevation sampled from a DEM raster.

`z_exterior` holds one elevation value per point in `base.exterior`.
`z_holes[k]` holds one elevation value per point in `base.holes[k]`.

`Geometry3D` objects are always derived from a `Geometry2D` + a DEM raster
via [`lift_to_3d`](@ref).  They are never constructed from native 3D vector
data (which does not exist in standard geospatial formats).
"""
struct Geometry3D
  base       :: Geometry2D
  z_exterior :: Vector{Float64}
  z_holes    :: Vector{Vector{Float64}}
end
