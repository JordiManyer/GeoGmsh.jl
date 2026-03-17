"""
Pyrenees 3D terrain mesh example.

Downloads:
  • NUTS-2 boundary for Navarre / Comunidad Foral de Navarra (ES22)
    from the Eurostat GISCO distribution service.
  • Copernicus GLO-30 DEM tiles (30 m resolution) for the same area
    from the public AWS S3 bucket (no authentication required).

Produces a terrain-following 3D surface mesh (.msh) whose nodes are lifted
to the actual elevation of the Pyrenean landscape.

DEM tiles used (1°×1° each, GLO-30 convention):
  Latitude band N41–N43, longitude columns W003–E000
  — covers roughly 41–44°N, 3°W–1°E, which encloses Navarre.

CRS workflow:
  The Copernicus tiles are natively in EPSG:4326 (geographic degrees).
  Both the vector boundary and the DEM are reprojected to UTM zone 30N
  (EPSG:25830, metres) so that mesh_size can be specified in metres and
  elevation values are physically consistent with horizontal coordinates.
"""

using GeoGmsh
import ArchGDAL
using Downloads
import GDAL_jll

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ---------------------------------------------------------------------------
# 1. Download NUTS-2 boundary
# ---------------------------------------------------------------------------

nuts_url  = "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/" *
            "NUTS_RG_01M_2024_4326_LEVL_2.geojson"
nuts_path = joinpath(data_dir, "NUTS_RG_01M_2024_4326_LEVL_2.geojson")

if !isfile(nuts_path)
  println("Downloading NUTS-2 boundaries (~16 MB)…")
  Downloads.download(nuts_url, nuts_path)
  println("  Saved: ", nuts_path)
end

# ---------------------------------------------------------------------------
# 2. Inspect
# ---------------------------------------------------------------------------

println("\nNavarre (ES22):")
comps = list_components(nuts_path; select = row -> row.NUTS_ID == "ES22")
println(comps)

# ---------------------------------------------------------------------------
# 3. Download Copernicus GLO-30 DEM tiles
# ---------------------------------------------------------------------------

# Each tile covers 1°×1° and is labelled by its south-west corner.
# Tiles below cover the full bounding box of Navarre (≈ 41.9–43.3°N, 2.2°W–0.7°E)
# with a one-tile margin on each side.
tiles = [
  ("N41_00", "W003_00"), ("N41_00", "W002_00"), ("N41_00", "W001_00"), ("N41_00", "E000_00"),
  ("N42_00", "W003_00"), ("N42_00", "W002_00"), ("N42_00", "W001_00"), ("N42_00", "E000_00"),
  ("N43_00", "W003_00"), ("N43_00", "W002_00"), ("N43_00", "W001_00"), ("N43_00", "E000_00"),
]

const DEM_BASE = "https://copernicus-dem-30m.s3.amazonaws.com"

println("\nDownloading DEM tiles…")
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

# ---------------------------------------------------------------------------
# 4. Mosaic the DEM tiles into a single VRT (EPSG:4326)
# ---------------------------------------------------------------------------

# gdalbuildvrt creates a seamless virtual raster from multiple GeoTIFFs.
# GDAL_jll is a transitive dependency of ArchGDAL.jl and provides the binary.
dem_vrt_4326 = joinpath(data_dir, "pyrenees_dem_4326.vrt")
GDAL_jll.gdalbuildvrt() do exe
  run(`$exe $dem_vrt_4326 $(tile_paths...)`)
end
println("\nMosaic (EPSG:4326): ", dem_vrt_4326)

# ---------------------------------------------------------------------------
# 5. Reproject DEM to UTM zone 30N (EPSG:25830, metres)
# ---------------------------------------------------------------------------

# The vector boundary will also be reprojected to EPSG:25830 by the pipeline.
# Both datasets must be in the same CRS before elevation sampling.
dem_tif_utm = joinpath(data_dir, "pyrenees_dem_25830.tif")
if !isfile(dem_tif_utm)
  println("Reprojecting DEM to EPSG:25830…")
  GDAL_jll.gdalwarp() do exe
    run(`$exe -t_srs EPSG:25830 -tr 30 30 -r bilinear
             -co COMPRESS=DEFLATE
             $dem_vrt_4326 $dem_tif_utm`)
  end
  println("  Saved: ", dem_tif_utm)
end

# ---------------------------------------------------------------------------
# 6. Terrain mesh
# ---------------------------------------------------------------------------

output = joinpath(data_dir, "pyrenees")

# target_crs = "EPSG:25830" reprojects the vector boundary to UTM 30N,
# matching the DEM.  mesh_size is now in metres (5 km element size).
geoms_to_msh_3d(
  nuts_path, dem_tif_utm, output;
  select         = row -> row.NUTS_ID == "ES22",
  target_crs     = "EPSG:25830",
  simplify_alg   = MinEdgeLength(tol = 5_000.0),
  mesh_size      = 5_000.0,    # 5 km characteristic element length
  nodata_fill    = 0.0,
  verbose        = true,
)
