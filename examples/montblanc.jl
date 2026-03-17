# # Mont Blanc — 3D terrain mesh
#
# This example produces two terrain-following meshes for the Mont Blanc massif
# using a user-defined bounding box and a Copernicus GLO-30 Digital Elevation
# Model:
#
# - **Surface mesh** (`geoms_to_msh_3d`): a single terrain-following surface
#   with z-coordinates sampled from the DEM.
# - **Volume mesh** (`geoms_to_msh_3d_volume`): the same surface extruded
#   downward by `depth` metres to form a solid tetrahedral block.
#
# **Features highlighted:**
# - Defining the domain as a bounding-box GeoJSON polygon (no administrative
#   boundary needed)
# - Downloading and mosaicking Copernicus GLO-30 DEM tiles with `GDAL_jll`
# - `geoms_to_msh_3d` and `geoms_to_msh_3d_volume`
# - Choosing a UTM CRS so that `mesh_size` and `depth` are in metres
#
# !!! note "Aspect ratio"
#     The Mont Blanc massif spans ~37 km × ~27 km horizontally but only
#     ~4,800 m vertically.  At true scale the mesh looks flat — apply vertical
#     exaggeration (e.g. 5×) in your visualiser (ParaView, QGIS 3D, Visit).
#
# | Mont Blanc (3D terrain) |
# |:-----------------------:|
# | ![Mont Blanc mesh](../assets/montblanc.png) |

using GeoGmsh
using Downloads
import GDAL_jll

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ## Bounding box
#
# We define the domain of interest as a simple rectangular polygon in
# EPSG:4326.  A small padding around the target area avoids edge effects
# in the mesh.
#
# Original bounding box: `6.785431, 45.785243, 7.044296, 45.956878`

lon_min, lat_min = 6.73, 45.75
lon_max, lat_max = 7.10, 45.99

bbox_path = joinpath(data_dir, "montblanc_bbox.geojson")
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
# The Copernicus GLO-30 DEM is distributed as 1°×1° GeoTIFF tiles labelled
# by their south-west corner.  Tiles are freely available from the public AWS
# S3 bucket — no authentication required.
#
# Our domain sits in the N45 latitude band (45–46°N) and spans longitude
# columns E006 and E007.

tiles = [("N45_00", "E006_00"), ("N45_00", "E007_00")]

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
#
# `gdalbuildvrt` stitches the individual tiles into a seamless virtual raster
# (VRT).  `gdalwarp` then reprojects it to UTM zone 32N (EPSG:32632) at 30 m
# resolution so that coordinates and elevations share the same unit (metres).

dem_vrt_4326 = joinpath(data_dir, "montblanc_dem_4326.vrt")
GDAL_jll.gdalbuildvrt_exe() do exe
  run(`$exe $dem_vrt_4326 $tile_paths`)
end
println("Mosaic (EPSG:4326): ", dem_vrt_4326)

dem_tif_utm = joinpath(data_dir, "montblanc_dem_32632.tif")
if !isfile(dem_tif_utm)
  println("Reprojecting DEM to EPSG:32632…")
  GDAL_jll.gdalwarp_exe() do exe
    run(`$exe -t_srs EPSG:32632 -tr 30 30 -r bilinear
             -co COMPRESS=DEFLATE
             $dem_vrt_4326 $dem_tif_utm`)
  end
  println("  Saved: ", dem_tif_utm)
end

# ## 3D terrain meshes
#
# `geoms_to_msh_3d` runs the standard 2D pipeline (reproject → simplify →
# ingest) and then lifts every mesh node's z-coordinate by bilinearly
# interpolating the DEM at its (x, y) position.
# `mesh_size = 500.0` gives ~500 m characteristic element length (in metres,
# consistent with the UTM CRS).
#
# `geoms_to_msh_3d_volume` does the same but also extrudes the surface
# downward by `depth` metres to produce a solid tetrahedral mesh.
# A depth of 1,000 m gives a thin but clearly visible pedestal relative to the
# ~4,800 m of vertical relief in the massif.

output = joinpath(data_dir, "montblanc")

geoms_to_msh_3d(
  bbox_path, dem_tif_utm, output;
  target_crs   = "EPSG:32632",
  simplify_alg = MinEdgeLength(tol = 500.0),
  mesh_size    = 500.0,
  nodata_fill  = 0.0,
  verbose      = true,
)

geoms_to_msh_3d_volume(
  bbox_path, dem_tif_utm, output * "_volume";
  target_crs   = "EPSG:32632",
  simplify_alg = MinEdgeLength(tol = 500.0),
  mesh_size    = 500.0,
  depth        = 1_000.0,
  nodata_fill  = 0.0,
  verbose      = true,
)
