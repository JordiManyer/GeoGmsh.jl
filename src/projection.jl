"""
Geometry rescaling utilities.

## Rescaling

    rescale(geoms, L)

Uniformly scales and translates `geoms` so that the largest bounding-box
dimension becomes exactly `L` and the minimum corner sits at the origin.

Reprojection is now handled upstream by `GeometryOps.reproject` before
geometries are ingested.
"""

# ============================================================================
# Public API
# ============================================================================

"""
    rescale(geoms, L) -> Vector{Geometry2D}

Uniformly scale and translate `geoms` so that the largest dimension of the
global bounding box equals `L` and the minimum corner is at the origin.
"""
function rescale(geoms::Vector{Geometry2D}, L::Real) :: Vector{Geometry2D}
  xmin, xmax, ymin, ymax = _global_bbox(geoms)
  scale = Float64(L) / max(xmax - xmin, ymax - ymin)
  f = pt -> ((pt[1] - xmin) * scale, (pt[2] - ymin) * scale)
  return [_project_geometry(g, f) for g in geoms]
end

# ============================================================================
# Internal helpers
# ============================================================================

function _project_contour(c::Contour, f)
  Contour(map(f, c.points), c.closed)
end

function _project_geometry(g::Geometry2D, f)
  Geometry2D(
    _project_contour(g.exterior, f),
    [_project_contour(h, f) for h in g.holes],
    g.name,
  )
end

function _global_bbox(geoms::Vector{Geometry2D})
  xmin = ymin =  Inf
  xmax = ymax = -Inf
  for g in geoms
    for pt in g.exterior.points
      xmin = min(xmin, pt[1]); xmax = max(xmax, pt[1])
      ymin = min(ymin, pt[2]); ymax = max(ymax, pt[2])
    end
    for h in g.holes, pt in h.points
      xmin = min(xmin, pt[1]); xmax = max(xmax, pt[1])
      ymin = min(ymin, pt[2]); ymax = max(ymax, pt[2])
    end
  end
  return xmin, xmax, ymin, ymax
end
