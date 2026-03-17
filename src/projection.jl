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

# ============================================================================
# Auto UTM helpers
# ============================================================================

# Return the centre of the bounding box in source CRS units.
# Two methods: one for DataFrames, one for raw GI geometries.
function _bbox_centre_lonlat(df::DataFrames.AbstractDataFrame)
  col  = first(GI.geometrycolumns(df))
  xmin = ymin =  Inf
  xmax = ymax = -Inf
  for g in df[!, col]
    isnothing(g) && continue
    ext = GI.extent(g)
    isnothing(ext) && continue
    xmin = min(xmin, ext.X[1]);  xmax = max(xmax, ext.X[2])
    ymin = min(ymin, ext.Y[1]);  ymax = max(ymax, ext.Y[2])
  end
  (isinf(xmin) || isinf(ymin)) && return 0.0, 0.0
  return (xmin + xmax) / 2.0, (ymin + ymax) / 2.0
end

function _bbox_centre_lonlat(geom)
  ext = GI.extent(geom)
  isnothing(ext) && return 0.0, 0.0
  return (ext.X[1] + ext.X[2]) / 2.0, (ext.Y[1] + ext.Y[2]) / 2.0
end

# Derive the UTM EPSG string from a geometry's bounding-box centre.
# Assumes the geometry is in geographic coordinates (lon/lat, EPSG:4326).
function _auto_utm_crs(geom) :: String
  lon, lat = _bbox_centre_lonlat(geom)
  zone = clamp(floor(Int, (lon + 180.0) / 6.0) + 1, 1, 60)
  ns   = lat >= 0.0 ? "6" : "7"
  return "EPSG:32" * ns * lpad(zone, 2, '0')
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
