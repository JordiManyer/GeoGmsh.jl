# # Everest — 3D terrain mesh
#
# This example produces terrain-following 3D surface and volume meshes for the
# Mount Everest region using a user-defined bounding box and Copernicus GLO-30
# DEM tiles.  The workflow is identical to the Mont Blanc example but covers a
# larger area (four DEM tiles) and uses UTM zone 45N.
#
# **Features highlighted:**
# - `prepare_dem` with a four-tile domain
# - `geoms_to_msh_3d` and `geoms_to_msh_3d_volume`
#
# !!! note "Aspect ratio"
#     The domain spans ~84 km × ~66 km horizontally; Everest is 8,849 m tall.
#     The true aspect ratio is roughly 1:9 — apply vertical exaggeration in
#     your visualiser to make the topography visible.
#
# | Everest (3D terrain) |
# |:--------------------:|
# | ![Everest mesh](../assets/everest.png) |

using GeoGmsh

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ## Bounding box
#
# Domain of interest with a small padding:
# original box: `86.551666, 27.712710, 87.301483, 28.214870`

lon_min, lat_min = 86.50, 27.66
lon_max, lat_max = 87.35, 28.26
bbox = (lon_min, lat_min, lon_max, lat_max)

bbox_path = joinpath(data_dir, "everest_bbox.geojson")
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

# ## Prepare DEM
#
# Four 1°×1° tiles cover the padded domain: latitude bands N27 and N28,
# longitude columns E086 and E087.  `prepare_dem` downloads, mosaics, and
# reprojects them to UTM zone 45N (EPSG:32645) in one call.

dem_tif_utm = prepare_dem(
  bbox, "EPSG:32645",
  joinpath(data_dir, "everest_dem_32645.tif");
  source   = :cop30,
  dest_dir = data_dir,
)

# ## 3D terrain meshes
#
# A depth of 1,000 m gives a thin but clearly visible pedestal relative to the
# ~6,000 m of vertical relief in the domain.

output = joinpath(data_dir, "everest")

geoms_to_msh_3d(
  bbox_path, dem_tif_utm, output;
  target_crs   = "EPSG:32645",
  simplify_alg = MinEdgeLength(tol = 500.0),
  mesh_size    = 500.0,
  nodata_fill  = 0.0,
  verbose      = true,
)

geoms_to_msh_3d_volume(
  bbox_path, dem_tif_utm, output * "_volume";
  target_crs   = "EPSG:32645",
  simplify_alg = MinEdgeLength(tol = 500.0),
  mesh_size    = 500.0,
  depth        = 1_000.0,
  nodata_fill  = 0.0,
  verbose      = true,
)
