"""
Coordinate projection and geometry rescaling utilities.

## Reprojection

    project_to_meters(geoms, source_crs; target = "EPSG:3857")

Uses the PROJ library (via Proj.jl) to reproject `geoms` from `source_crs`
to `target`.  `source_crs` is the raw WKT string returned by `read_shapefile`
(or any other CRS string PROJ understands).  `target` may be a CRS string
(e.g. `"EPSG:3857"`) or a pre-built `Proj.Transformation`.

## Rescaling

    rescale(geoms, L)

Uniformly scales and translates `geoms` so that the largest bounding-box
dimension becomes exactly `L` and the minimum corner sits at the origin.
"""

import Proj

# ============================================================================
# Public API
# ============================================================================

"""
    project_to_meters(geoms, source_crs; target = "EPSG:3857") -> Vector{ShapeGeometry}

Reproject `geoms` using the PROJ library.

- `source_crs` — WKT (or any PROJ-recognised CRS string) for the input data,
  as returned by `read_shapefile`.  Pass `nothing` to assume EPSG:4326
  (geographic degrees) with a warning.
- `target` — destination CRS string (default `"EPSG:3857"`) or a pre-built
  `Proj.Transformation`.
"""
function project_to_meters(
  geoms      :: Vector{ShapeGeometry},
  source_crs :: Union{String,Nothing};
  target     :: Union{String,Proj.Transformation} = "EPSG:3857",
) :: Vector{ShapeGeometry}
  if isnothing(source_crs)
    @warn "No CRS information found; assuming EPSG:4326 (geographic degrees)."
    source_crs = "EPSG:4326"
  end

  trans = target isa Proj.Transformation ? target :
          Proj.Transformation(source_crs, target; always_xy = true)

  project = function(pt::NTuple{2,Float64})
    r = trans(pt[1], pt[2])
    (Float64(r[1]), Float64(r[2]))
  end

  return [_project_geometry(g, project) for g in geoms]
end

"""
    rescale(geoms, L) -> Vector{ShapeGeometry}

Uniformly scale and translate `geoms` so that the largest dimension of the
global bounding box equals `L` and the minimum corner is at the origin.
"""
function rescale(geoms::Vector{ShapeGeometry}, L::Real) :: Vector{ShapeGeometry}
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

function _project_geometry(g::ShapeGeometry, f)
  ShapeGeometry(
    _project_contour(g.exterior, f),
    [_project_contour(h, f) for h in g.holes],
    g.name,
  )
end

function _global_bbox(geoms::Vector{ShapeGeometry})
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
