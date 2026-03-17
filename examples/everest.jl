# # Everest — 3D terrain mesh
#
# This example produces terrain-following 3D surface and volume meshes for the
# Mount Everest region using a user-defined bounding box and Copernicus GLO-30
# DEM tiles.  The workflow is identical to the Mont Blanc example but covers a
# larger area (four DEM tiles) and uses UTM zone 45N.
#
# **Features highlighted:**
# - 3D surface and volume meshing with `geoms_to_msh_3d` and
#   `geoms_to_msh_3d_volume` over a larger multi-tile domain
# - Choosing the correct UTM zone for the area of interest (zone 45N for Nepal)
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
using Downloads
import GDAL_jll

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ## Bounding box
#
# Domain of interest with a small padding:
# original box: `86.551666, 27.712710, 87.301483, 28.214870`

lon_min, lat_min = 86.50, 27.66
lon_max, lat_max = 87.35, 28.26

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

# ## Download DEM tiles
#
# Four 1°×1° tiles cover the padded domain:
# latitude bands N27 and N28, longitude columns E086 and E087.

tiles = [
  ("N27_00", "E086_00"), ("N27_00", "E087_00"),
  ("N28_00", "E086_00"), ("N28_00", "E087_00"),
]

const DEM_BASE = "https://copernicus-dem-30m.s3.amazonaws.com"

println("Downloading DEM tiles…")
tile_paths = String[]
for (lat, lon) in tiles
  name = "Copernicus_DSM_COG_10_$(lat)_$(lon)_DEM"
  url  = "$DEM_BASE/$name/$name.tif"
  dest = joinpath(data_dir, "$name.tif")
  if isfile(dest)
    print("  $name  (cached)\n")
  else
    print("  $name  downloading… ")
    try
      Downloads.download(url, dest)
      println("ok")
    catch e
      println("FAILED: $e")
      continue
    end
  end
  push!(tile_paths, dest)
end
println("$(length(tile_paths)) tile(s) ready.")

# ## Mosaic and reproject

dem_vrt_4326 = joinpath(data_dir, "everest_dem_4326.vrt")
GDAL_jll.gdalbuildvrt_exe() do exe
  run(`$exe $dem_vrt_4326 $tile_paths`)
end
println("Mosaic (EPSG:4326): ", dem_vrt_4326)

dem_tif_utm = joinpath(data_dir, "everest_dem_32645.tif")
if !isfile(dem_tif_utm)
  println("Reprojecting DEM to EPSG:32645…")
  GDAL_jll.gdalwarp_exe() do exe
    run(`$exe -t_srs EPSG:32645 -tr 30 30 -r bilinear
             -co COMPRESS=DEFLATE
             $dem_vrt_4326 $dem_tif_utm`)
  end
  println("  Saved: ", dem_tif_utm)
end

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
