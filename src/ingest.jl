"""
    ingest(geom) -> Vector{Geometry2D}

Convert any GeoInterface-compatible geometry to our internal representation.

Accepted inputs:
- Any `GI.PolygonTrait` geometry (single polygon)
- Any `GI.MultiPolygonTrait` geometry (split into one `Geometry2D` each)
- Any `GI.FeatureTrait` (geometry is extracted; properties are dropped)
- Any `GI.FeatureCollectionTrait` (all features ingested and concatenated)
- A `DataFrame` with a geometry column (as returned by `read_geodata`)

Ring orientation is normalised: exterior rings → CCW, hole rings → CW.
Coordinates are always converted to `Float64`.
Closing points (repeated first point) are removed automatically.
"""
function ingest(geom) :: Vector{Geometry2D}
  # GeoInterface uses isfeature/isfeaturecollection rather than geomtrait
  # for Feature and FeatureCollection types.
  if GI.isfeaturecollection(geom)
    return _ingest_fc(geom)
  elseif GI.isfeature(geom)
    return _ingest_feat(geom)
  else
    return _ingest(GI.geomtrait(geom), geom)
  end
end

function ingest(df::DataFrames.AbstractDataFrame) :: Vector{Geometry2D}
  col = first(GI.geometrycolumns(df))
  result = Geometry2D[]
  for row in eachrow(df)
    g = row[col]
    isnothing(g) && continue
    append!(result, ingest(g))
  end
  return result
end

# ---------------------------------------------------------------------------
# Trait dispatch
# ---------------------------------------------------------------------------

function _ingest(::GI.PolygonTrait, geom)
  Geometry2D[_polygon_to_shape(geom)]
end

function _ingest(::GI.MultiPolygonTrait, geom)
  [_polygon_to_shape(GI.getgeom(geom, i)) for i in 1:GI.ngeom(geom)]
end

function _ingest_feat(feat)
  g = GI.geometry(feat)
  isnothing(g) ? Geometry2D[] : ingest(g)
end

function _ingest_fc(fc)
  result = Geometry2D[]
  for i in 1:GI.nfeature(fc)
    append!(result, ingest(GI.getfeature(fc, i)))
  end
  return result
end

function _ingest(trait, geom)
  @warn "ingest: unsupported geometry type $(typeof(trait)); skipping."
  return Geometry2D[]
end

# ---------------------------------------------------------------------------
# Ring → Contour conversion
# ---------------------------------------------------------------------------

# Convert any GI ring (LinearRingTrait or LineStringTrait) to a
# normalised Contour with the requested orientation (:ccw or :cw).
function _polygon_to_shape(poly) :: Geometry2D
  ext_pts  = _pts_from_ring(GI.getexterior(poly))
  hole_pts = [_pts_from_ring(GI.gethole(poly, i)) for i in 1:GI.nhole(poly)]

  exterior = Contour(_orient(ext_pts, :ccw), true)
  holes    = [Contour(_orient(h, :cw), true) for h in hole_pts
              if length(h) >= 3]

  return Geometry2D(exterior, holes)
end

# Extract NTuple{2,Float64} points from any GI ring, stripping the
# explicit closing point if the ring repeats its first point.
function _pts_from_ring(ring) :: Vector{NTuple{2,Float64}}
  n   = GI.npoint(ring)
  pts = NTuple{2,Float64}[
    (Float64(GI.x(GI.getpoint(ring, i))), Float64(GI.y(GI.getpoint(ring, i))))
    for i in 1:n
  ]
  # GI rings (and ArchGDAL LineStrings used as rings) often repeat first point.
  if length(pts) > 1 && pts[1] == pts[end]
    pop!(pts)
  end
  return pts
end

# Reverse pts if they don't match the requested orientation.
# :ccw → positive signed area; :cw → negative signed area.
function _orient(pts::Vector{NTuple{2,Float64}}, orientation::Symbol)
  area = _ring_signed_area(pts)
  if orientation == :ccw && area < 0
    reverse!(pts)
  elseif orientation == :cw && area > 0
    reverse!(pts)
  end
  return pts
end

# Shoelace formula on raw points. Positive = CCW, Negative = CW.
function _ring_signed_area(pts::Vector{NTuple{2,Float64}})
  n = length(pts)
  n < 3 && return 0.0
  a = 0.0
  for i in 1:n
    j = mod1(i + 1, n)
    a += pts[i][1] * pts[j][2] - pts[j][1] * pts[i][2]
  end
  return a / 2
end
