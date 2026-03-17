"""
    read_geodata(path; layer = nothing, select = nothing) -> DataFrame

Read any geospatial file supported by GeoDataFrames (Shapefile, GeoJSON,
GeoPackage, GeoParquet, GeoArrow, FlatGeobuf, and any GDAL-supported format).

Returns a standard `DataFrame` with a `:geometry` column of
GeoInterface-compatible geometries. CRS metadata is accessible via
`GeoInterface.crs(df)`.

# Keyword arguments
- `layer`  — layer index (`Int`, 0-based) or name (`String`) for multi-layer
             formats such as GeoPackage.  Defaults to the first layer.
- `select` — a predicate `row -> Bool` to filter rows after reading.
"""
function read_geodata(path::AbstractString; layer = nothing, select = nothing)
  df = isnothing(layer) ? GeoDataFrames.read(path) :
                          GeoDataFrames.read(path; layer)
  df = _expand_rings(df)
  isnothing(select) ? df : filter(select, df)
end

"""
    list_components(path; layer = nothing) -> DataFrame

Read the geospatial file at `path` and print a summary table that includes
all attribute columns plus per-geometry statistics (`:n_pts`, `:area`,
`:xmin`, `:xmax`, `:ymin`, `:ymax`). Returns the augmented DataFrame.

Useful for exploring a file before deciding which features to select.
For multi-layer formats, pass `layer` to inspect a specific layer.
"""
function list_components(path::AbstractString; layer = nothing)
  df = read_geodata(path; layer)   # :ring, :n_pts, :area, :xmin/:xmax/:ymin/:ymax already present
  isempty(df) && (println("(empty)"); return df)
  show(df; allcols = true)
  println()
  return df
end

"""
    read_shapefile(path; layer = nothing, select = nothing) -> DataFrame

Thin backward-compatible wrapper around [`read_geodata`](@ref).
"""
read_shapefile(path::AbstractString; kwargs...) = read_geodata(path; kwargs...)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Expand MultiPolygon rows into one row per polygon component, sorted by area
# descending.  Adds :ring (1 = largest), :n_pts, :area, :xmin/:xmax/:ymin/:ymax.
# Polygon rows get ring = 1.  Rows with missing geometry get ring = 0 and NaN stats.
# Column-level metadata (e.g. CRS) on the geometry column is preserved.
function _expand_rings(df::DataFrames.AbstractDataFrame)
  isempty(df) && return df
  geom_col = first(GI.geometrycolumns(df))

  # transform adds vector-valued columns (one entry per component); flatten expands them.
  out      = transform(df, geom_col => ByRow(_geom_components) => AsTable)
  crs_meta = DataFrames.colmetadata(out, geom_col)   # capture before drop
  select!(out, Not(geom_col))
  rename!(out, :_geom => geom_col)
  for (k, v) in crs_meta                             # restore CRS on the new column
    DataFrames.colmetadata!(out, geom_col, k, v; style = :note)
  end
  return flatten(out, [geom_col, :ring, :n_pts, :area, :xmin, :xmax, :ymin, :ymax])
end

# Statistics for a single polygon's exterior ring.
function _exterior_stats(poly)
  pts = _pts_from_ring(GI.getexterior(poly))
  isempty(pts) && return (n_pts=0, area=NaN, xmin=NaN, xmax=NaN, ymin=NaN, ymax=NaN)
  xs = [p[1] for p in pts]
  ys = [p[2] for p in pts]
  (n_pts = length(pts),
   area  = abs(_ring_signed_area(pts)),
   xmin  = minimum(xs), xmax = maximum(xs),
   ymin  = minimum(ys), ymax = maximum(ys))
end

# Per-geometry exploder: returns a NamedTuple of aligned vectors, one entry per
# polygon component.  Used by _expand_rings via ByRow + AsTable + flatten.
function _geom_components(g)
  if isnothing(g)
    return (_geom=[g], ring=[0], n_pts=[0], area=[NaN],
            xmin=[NaN], xmax=[NaN], ymin=[NaN], ymax=[NaN])
  end
  if GI.geomtrait(g) isa GI.MultiPolygonTrait
    polys = [GI.getgeom(g, j) for j in 1:GI.ngeom(g)]
    stats = [_exterior_stats(p) for p in polys]
    order = sortperm([s.area for s in stats], rev = true)
    return (_geom = [polys[i] for i in order],
            ring  = collect(1:length(polys)),
            n_pts = [stats[i].n_pts for i in order],
            area  = [stats[i].area  for i in order],
            xmin  = [stats[i].xmin  for i in order],
            xmax  = [stats[i].xmax  for i in order],
            ymin  = [stats[i].ymin  for i in order],
            ymax  = [stats[i].ymax  for i in order])
  else
    s = _exterior_stats(g)
    return (_geom=[g], ring=[1], n_pts=[s.n_pts], area=[s.area],
            xmin=[s.xmin], xmax=[s.xmax], ymin=[s.ymin], ymax=[s.ymax])
  end
end

# Extract exterior ring points from any GI polygon or multipolygon.
# For multipolygons, uses the largest polygon's exterior.
function _all_exterior_pts(geom)
  t = GI.geomtrait(geom)
  _all_exterior_pts(t, geom)
end

function _all_exterior_pts(::GI.PolygonTrait, geom)
  _pts_from_ring(GI.getexterior(geom))
end

function _all_exterior_pts(::GI.MultiPolygonTrait, geom)
  pts = NTuple{2,Float64}[]
  for i in 1:GI.ngeom(geom)
    append!(pts, _pts_from_ring(GI.getexterior(GI.getgeom(geom, i))))
  end
  return pts
end

function _all_exterior_pts(::GI.FeatureTrait, feat)
  g = GI.geometry(feat)
  isnothing(g) ? NTuple{2,Float64}[] : _all_exterior_pts(g)
end

function _all_exterior_pts(_, _)
  NTuple{2,Float64}[]
end

# Extract the WKT string from whatever GI.crs() returns so that
# Proj.Transformation can consume it.
function _crs_to_wkt(crs) :: Union{String, Nothing}
  isnothing(crs) && return nothing
  crs isa AbstractString && return String(crs)
  # GeoFormatTypes objects expose the underlying value via .val
  return string(crs.val)
end
