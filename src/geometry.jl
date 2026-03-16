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
    ShapeGeometry

A single polygon: an exterior ring plus zero or more hole rings.
All rings are stored as `Contour` with `closed = true`.
Exterior rings are oriented CCW (positive signed area);
hole rings are oriented CW (negative signed area).
"""
struct ShapeGeometry
  exterior::Contour
  holes::Vector{Contour}
  name::String
end

ShapeGeometry(exterior::Contour, holes::Vector{Contour}) =
  ShapeGeometry(exterior, holes, "")

# ---------------------------------------------------------------------------
# GeoInterface traits
# ---------------------------------------------------------------------------

import GeoInterface as GI

GI.isgeometry(::Type{ShapeGeometry}) = true
GI.geomtrait(::ShapeGeometry)        = GI.PolygonTrait()
GI.nhole(::GI.PolygonTrait, g::ShapeGeometry)        = length(g.holes)
GI.getexterior(::GI.PolygonTrait, g::ShapeGeometry)  = g.exterior
GI.gethole(::GI.PolygonTrait, g::ShapeGeometry, i)   = g.holes[i]

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
function npoints(g::ShapeGeometry)
  n = npoints(g.exterior)
  for h in g.holes
    n += npoints(h)
  end
  return n
end
