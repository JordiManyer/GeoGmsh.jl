"""
Terrain-aware 2+1D meshing support.

Provides tools to read Digital Elevation Model (DEM) rasters and lift 2D
polygon geometries into 3D by sampling elevation at every boundary point.

## Workflow

```julia
dem   = read_dem("srtm_tile.tif")
geoms = ingest(read_geodata("region.shp"))   # Vector{Geometry2D}
g3d   = lift_to_3d(geoms, dem)              # Vector{Geometry3D}
write_geo(g3d, "output"; mesh_size = 1.0)
```

## CRS note

The DEM and the 2D geometries must be in the same CRS before sampling.
Use `GeometryOps.reproject` (via `target_crs` in `geoms_to_geo_3d`) to
reproject the vector data, and ensure the DEM is in a matching projection.
"""

import ArchGDAL

# ============================================================================
# DEMRaster
# ============================================================================

"""
    DEMRaster

An in-memory Digital Elevation Model raster.

Fields:
- `data`      — elevation matrix, shape `(ncols, nrows)` (ArchGDAL convention).
                 `data[col, row]` where col=1 is west, row=1 is north.
- `transform` — GDAL 6-element geotransform `[x_origin, x_pixel, x_rot,
                 y_origin, y_rot, y_pixel]`.  For north-up rasters `y_pixel < 0`.
- `crs_wkt`   — WKT CRS string, or `nothing` if unavailable.
- `nodata`    — nodata sentinel value, or `nothing` if not set.
"""
struct DEMRaster
  data      :: Matrix{Float64}
  transform :: Vector{Float64}
  crs_wkt   :: Union{String,Nothing}
  nodata    :: Union{Float64,Nothing}
end

# ============================================================================
# Public API
# ============================================================================

"""
    read_dem(path) -> DEMRaster

Read a Digital Elevation Model raster file (GeoTIFF, SRTM .hgt, NetCDF, …)
using ArchGDAL. Only the first band is read.

Returns a [`DEMRaster`](@ref) with the elevation data, geotransform, CRS, and
nodata value.
"""
function read_dem(path::AbstractString) :: DEMRaster
  ArchGDAL.read(path) do dataset
    band      = ArchGDAL.getband(dataset, 1)
    data      = Float64.(ArchGDAL.read(band))   # (ncols, nrows)
    transform = Float64.(ArchGDAL.getgeotransform(dataset))
    crs_wkt   = let p = ArchGDAL.getproj(dataset); isempty(p) ? nothing : p end
    nodata    = ArchGDAL.getnodatavalue(band)
    DEMRaster(data, transform, crs_wkt, nodata)
  end
end

"""
    sample_elevation(pts, dem; nodata_fill = 0.0) -> Vector{Float64}

Bilinearly interpolate elevation values at `pts` (a `Vector{NTuple{2,Float64}}`)
from `dem`. Points that fall outside the raster extent or coincide with nodata
cells are filled with `nodata_fill` (default `0.0`).
"""
function sample_elevation(
  pts         :: Vector{NTuple{2,Float64}},
  dem         :: DEMRaster;
  nodata_fill :: Float64 = 0.0,
) :: Vector{Float64}
  gt    = dem.transform
  data  = dem.data
  ncols = size(data, 1)
  nrows = size(data, 2)
  nd    = dem.nodata

  result = Vector{Float64}(undef, length(pts))
  for (i, (x, y)) in enumerate(pts)
    # World coords → fractional pixel index (0-indexed, GDAL convention).
    col_f = (x - gt[1]) / gt[2]
    row_f = (y - gt[4]) / gt[6]

    # Integer pixel corners (0-indexed).
    c0 = floor(Int, col_f)
    r0 = floor(Int, row_f)

    # Out-of-bounds → fill.
    if c0 < 0 || c0 + 1 >= ncols || r0 < 0 || r0 + 1 >= nrows
      result[i] = nodata_fill
      continue
    end

    # Sample four neighbours (1-indexed Julia, data is [col, row]).
    z00 = data[c0 + 1, r0 + 1]
    z10 = data[c0 + 2, r0 + 1]
    z01 = data[c0 + 1, r0 + 2]
    z11 = data[c0 + 2, r0 + 2]

    # Nodata check.
    if !isnothing(nd) && any(z -> z == nd, (z00, z10, z01, z11))
      result[i] = nodata_fill
      continue
    end

    # Bilinear weights.
    tc = col_f - c0
    tr = row_f - r0
    result[i] = (1-tc)*(1-tr)*z00 + tc*(1-tr)*z10 +
                (1-tc)*tr    *z01 + tc*tr    *z11
  end
  return result
end

"""
    lift_to_3d(geom,  dem; nodata_fill = 0.0) -> Geometry3D
    lift_to_3d(geoms, dem; nodata_fill = 0.0) -> Vector{Geometry3D}

Sample elevation from `dem` at every boundary point of `geom` (or each element
of `geoms`) and return the corresponding [`Geometry3D`](@ref) object(s).

The DEM and the 2D geometries must be in the same CRS.
"""
function lift_to_3d(g::Geometry2D, dem::DEMRaster; nodata_fill::Float64 = 0.0) :: Geometry3D
  z_ext   = sample_elevation(g.exterior.points, dem; nodata_fill)
  z_holes = [sample_elevation(h.points, dem; nodata_fill) for h in g.holes]
  Geometry3D(g, z_ext, z_holes)
end

function lift_to_3d(
  geoms       :: Vector{Geometry2D},
  dem         :: DEMRaster;
  nodata_fill :: Float64 = 0.0,
) :: Vector{Geometry3D}
  [lift_to_3d(g, dem; nodata_fill) for g in geoms]
end
