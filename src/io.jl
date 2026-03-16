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
  df = read_geodata(path; layer)
  isempty(df) && (println("(empty)"); return df)

  # Compute per-geometry statistics from the raw GI rings.
  n_pts = Int[]
  areas = Float64[]
  xmins = Float64[]
  xmaxs = Float64[]
  ymins = Float64[]
  ymaxs = Float64[]

  geom_col = first(GI.geometrycolumns(df))
  for row in eachrow(df)
    g = row[geom_col]
    if isnothing(g)
      push!(n_pts, 0); push!(areas, NaN)
      push!(xmins, NaN); push!(xmaxs, NaN)
      push!(ymins, NaN); push!(ymaxs, NaN)
      continue
    end
    pts = _all_exterior_pts(g)
    push!(n_pts, length(pts))
    push!(areas, isempty(pts) ? NaN : abs(_ring_signed_area(pts)))
    if isempty(pts)
      push!(xmins, NaN); push!(xmaxs, NaN)
      push!(ymins, NaN); push!(ymaxs, NaN)
    else
      xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
      push!(xmins, minimum(xs)); push!(xmaxs, maximum(xs))
      push!(ymins, minimum(ys)); push!(ymaxs, maximum(ys))
    end
  end

  df2 = copy(df)
  df2[!, :n_pts] = n_pts
  df2[!, :area]  = areas
  df2[!, :xmin]  = xmins
  df2[!, :xmax]  = xmaxs
  df2[!, :ymin]  = ymins
  df2[!, :ymax]  = ymaxs

  show(df2; allcols = true)
  println()
  return df2
end

"""
    read_shapefile(path; layer = nothing, select = nothing) -> DataFrame

Thin backward-compatible wrapper around [`read_geodata`](@ref).
"""
read_shapefile(path::AbstractString; kwargs...) = read_geodata(path; kwargs...)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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
