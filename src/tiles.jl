"""
DEM tile download and preparation utilities.

Provides:
- `download_cop30_tiles` — fetch Copernicus GLO-30 tiles for a bounding box.
- `prepare_dem`          — mosaic + reproject tile rasters (generic), with a
                           convenience overload that also downloads first.
"""

import GDAL_jll
using Downloads

# ============================================================================
# Copernicus GLO-30
# ============================================================================

const _COP30_BASE = "https://copernicus-dem-30m.s3.amazonaws.com"

function _cop30_tile_name(lat::Int, lon::Int) :: String
  lat_str = lat >= 0 ? "N$(lpad(lat, 2, '0'))_00" : "S$(lpad(-lat, 2, '0'))_00"
  lon_str = lon >= 0 ? "E$(lpad(lon, 3, '0'))_00" : "W$(lpad(-lon, 3, '0'))_00"
  return "Copernicus_DSM_COG_10_$(lat_str)_$(lon_str)_DEM"
end

"""
    download_cop30_tiles(bbox, dest_dir; verbose=true) -> Vector{String}

Download Copernicus GLO-30 DEM tiles covering `bbox = (lon_min, lat_min,
lon_max, lat_max)` from the public AWS S3 bucket into `dest_dir`.
Tiles already present on disk are skipped (cached).

Returns the paths of all successfully downloaded (or cached) tiles.
"""
function download_cop30_tiles(
  bbox     :: NTuple{4},
  dest_dir :: AbstractString;
  verbose  :: Bool = true,
) :: Vector{String}
  lon_min, lat_min, lon_max, lat_max = Float64.(bbox)
  mkpath(dest_dir)

  lats = floor(Int, lat_min) : floor(Int, lat_max)
  lons = floor(Int, lon_min) : floor(Int, lon_max)

  verbose && println("Downloading Copernicus GLO-30 tiles…")
  tile_paths = String[]
  for lat in lats, lon in lons
    name = _cop30_tile_name(lat, lon)
    url  = "$_COP30_BASE/$name/$name.tif"
    dest = joinpath(dest_dir, "$name.tif")
    if isfile(dest)
      verbose && println("  $name  (cached)")
    else
      verbose && print("  $name  downloading… ")
      try
        Downloads.download(url, dest)
        verbose && println("ok")
      catch e
        verbose && println("FAILED: $e")
        isfile(dest) && rm(dest)
        continue
      end
    end
    push!(tile_paths, dest)
  end
  verbose && println("$(length(tile_paths)) tile(s) ready.")
  return tile_paths
end

# ============================================================================
# prepare_dem
# ============================================================================

"""
    prepare_dem(tile_paths, target_crs, output_path; resolution=30.0, verbose=true)
    prepare_dem(bbox, target_crs, output_path; source=:cop30, dest_dir=..., ...)

Mosaic raster tiles and reproject to `target_crs`, writing a GeoTIFF to
`output_path`.  If `output_path` already exists it is returned immediately
(no work done — safe to call on every run).

The first form is **generic**: `tile_paths` can be any list of GDAL-readable
raster files (GeoTIFF, VRT, HGT, NetCDF, …).

The second form accepts a `bbox = (lon_min, lat_min, lon_max, lat_max)` and
downloads tiles automatically using the given `source` before mosaicking.

# Keyword arguments (both forms)
- `resolution` — output pixel size in `target_crs` units (default `30.0`; metres
                 for a projected UTM CRS, which matches the GLO-30 native grid).
- `verbose`    — print progress (default `true`).

# Additional keyword arguments (bbox form)
- `source`   — tile provider symbol; currently `:cop30` (Copernicus GLO-30,
               freely available from AWS S3).
- `dest_dir` — directory for downloaded tiles (default: directory of
               `output_path`, or the current working directory if `output_path`
               has no directory component).
"""
function prepare_dem(
  tile_paths  :: Vector{<:AbstractString},
  target_crs  :: AbstractString,
  output_path :: AbstractString;
  resolution  :: Real = 30.0,
  verbose     :: Bool = true,
) :: String
  if isfile(output_path)
    verbose && println("DEM cached: $output_path")
    return output_path
  end

  vrt = tempname() * ".vrt"
  try
    verbose && println("Mosaicking $(length(tile_paths)) tile(s)…")
    GDAL_jll.gdalbuildvrt_exe() do exe
      run(`$exe $vrt $tile_paths`)
    end

    res = Float64(resolution)
    verbose && println("Reprojecting → $target_crs at $(res) units/pixel…")
    GDAL_jll.gdalwarp_exe() do exe
      run(`$exe -t_srs $target_crs -tr $res $res -r bilinear
               -co COMPRESS=DEFLATE $vrt $output_path`)
    end
  finally
    isfile(vrt) && rm(vrt)
  end
  verbose && println("  Saved: $output_path")
  return output_path
end

function prepare_dem(
  bbox        :: NTuple{4},
  target_crs  :: AbstractString,
  output_path :: AbstractString;
  source      :: Symbol         = :cop30,
  dest_dir    :: AbstractString = let d = dirname(output_path); isempty(d) ? pwd() : d end,
  resolution  :: Real           = 30.0,
  verbose     :: Bool           = true,
) :: String
  if isfile(output_path)
    verbose && println("DEM cached: $output_path")
    return output_path
  end
  if source === :cop30
    tile_paths = download_cop30_tiles(bbox, dest_dir; verbose)
  else
    error("Unknown source :$source. Currently supported: :cop30")
  end
  prepare_dem(tile_paths, target_crs, output_path; resolution, verbose)
end
