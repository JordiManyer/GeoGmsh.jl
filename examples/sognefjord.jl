# # Sognefjord — 3D coastal/ocean-floor mesh
#
# This example produces a bidirectional tetrahedral volume mesh for the
# Sognefjord fjord system in western Norway, using a GEBCO bathymetry/
# topography DEM downloaded from the OpenTopography global API.
#
# **Features highlighted:**
# - `download_opentopo_dem`: one-call DEM retrieval from OpenTopography
# - `prepare_dem`: reproject the raw GeoTIFF to a projected UTM CRS
# - `geoms_to_msh_3d_volume` with `depth = (d_below, d_above)`:
#   bidirectional extrusion that captures both the water column and the
#   terrain above sea level in a single tetrahedral mesh
#
# **Domain:** roughly 4 °E – 6 °E, 60 °N – 62 °N (Sognefjord inlet area).
#
# !!! note "API key required"
#     `download_opentopo_dem` needs a free OpenTopography API key.
#     Register at https://opentopography.org/developers and either set
#     `ENV["OPENTOPO_API_KEY"]` or pass `api_key` explicitly.
#
# !!! note "GEBCO resolution"
#     GEBCO is a global ~500 m product.  For a mesh of the fjord system this
#     is sufficient; for finer coastal detail use `:cop30` (30 m) over land.

using GeoGmsh

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ## Bounding box
#
# The Sognefjord stretches roughly from 4 °E to 7 °E between 60 °N and 62 °N.
# We use a slightly smaller box centred on the fjord mouth to keep file sizes
# modest for the example.

lon_min, lat_min = 4.0, 60.0
lon_max, lat_max = 6.0, 62.0
bbox = (lon_min, lat_min, lon_max, lat_max)

bbox_path = joinpath(data_dir, "sognefjord_bbox.geojson")
open(bbox_path, "w") do io
  write(io, """
{
  "type": "FeatureCollection",
  "features": [{
    "type": "Feature",
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[$(lon_min),$(lat_min)],[$(lon_max),$(lat_min)],
                        [$(lon_max),$(lat_max)],[$(lon_min),$(lat_max)],
                        [$(lon_min),$(lat_min)]]]
    },
    "properties": {}
  }]
}
""")
end

# ## Download and prepare DEM
#
# `download_opentopo_dem` fetches the GEBCO Ice Surface Topography product
# (`:gebco`) from OpenTopography as a WGS-84 GeoTIFF.  The file is cached on
# disk — repeated runs skip the download.
#
# `prepare_dem` reprojects the raw tile to UTM zone 32N (EPSG:32632) at
# 500 m resolution, matching the native GEBCO grid spacing.

raw_tif = download_opentopo_dem(
  bbox, joinpath(data_dir, "sognefjord_gebco_raw.tif");
  demtype = :gebco,
  # api_key = "...",   # or set ENV["OPENTOPO_API_KEY"]
)

dem_tif = prepare_dem(
  [raw_tif], "EPSG:32632",
  joinpath(data_dir, "sognefjord_dem_32632.tif");
  resolution = 500.0,
)

# ## Bidirectional volume mesh
#
# `depth = (1_500.0, 500.0)` extrudes:
# - **1 500 m downward** — captures the deep fjord trench (max ~1 300 m)
# - **500 m upward**    — captures the coastal mountain shoulders
#
# The terrain surface becomes the "Interface" physical group in the mesh,
# separating the two sub-volumes ("Volume_Below" and "Volume_Above").
# "Bottom", "Top", and "Sides" are the outer boundary faces.

geoms_to_geo_3d(
  bbox_path, dem_tif,
  joinpath(data_dir, "sognefjord_surface");
  target_crs   = "EPSG:32632",
  simplify_alg = MinEdgeLength(tol = 500.0),
  mesh_size    = 500.0,
  nodata_fill  = 0.0,
  verbose      = true,
)

geoms_to_msh_3d_volume(
  bbox_path, dem_tif,
  joinpath(data_dir, "sognefjord_volume");
  target_crs   = "EPSG:32632",
  simplify_alg = MinEdgeLength(tol = 500.0),
  mesh_size    = 500.0,
  depth        = (500.0, 100.0),
  nodata_fill  = 0.0,
  verbose      = true,
)
