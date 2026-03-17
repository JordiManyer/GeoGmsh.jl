"""
    GeoGmsh

Convert geospatial data to Gmsh geometry (`.geo`) and mesh (`.msh`) files.
Accepts any GeoInterface-compatible source: Shapefiles, GeoJSON, GeoPackage,
GeoParquet, NaturalEarth data, or raw GeoInterface geometries.

## Typical workflow

```julia
using GeoGmsh

# From Natural Earth (no files needed)
using NaturalEarth
geoms_to_geo(naturalearth("admin_0_countries", 110), "countries";
  # target_crs defaults to :auto_utm — UTM zone detected from geometry centroid
  simplify_alg = MinEdgeLength(tol = 5_000.0),
  mesh_size    = 2.0,
)

# From any file (Shapefile, GeoJSON, GeoPackage, …)
df = read_geodata("regions.shp")
geoms_to_geo(df, "output"; mesh_size = 2.0)   # auto UTM by default

# 3D terrain mesh
geoms_to_msh_3d("region.shp", "srtm.tif", "terrain";
  target_crs = "EPSG:32632",
  mesh_size  = 500.0,
)
```

## Pipeline (2D)

`geoms_to_geo` / `geoms_to_msh` run these steps:

1. `GeometryOps.reproject`                — reproject coordinates (Proj.jl)
2. `GeometryOps.simplify(MinEdgeLength…)` — remove short edges
3. `GeometryOps.segmentize`               — split long edges
4. [`ingest`](@ref)                       — convert to Gmsh-ready `Geometry2D`
5. [`filter_components`](@ref)            — drop degenerate rings
6. [`rescale`](@ref)                      — normalise into an `L × L` bounding box
7. [`write_geo`](@ref) or [`generate_mesh`](@ref) — write output

## Pipeline (3D surface)

`geoms_to_geo_3d` / `geoms_to_msh_3d` add after step 6:

7. [`read_dem`](@ref)    — read DEM raster (GeoTIFF, SRTM, NetCDF, …)
8. [`lift_to_3d`](@ref)  — sample boundary point elevations
9. [`write_geo`](@ref) or [`generate_mesh`](@ref) with [`Geometry3D`](@ref)

## Pipeline (3D volume)

`geoms_to_msh_3d_volume` produces a **tetrahedral** solid mesh by extruding
the terrain surface downward by `depth` (in CRS units):

7. [`read_dem`](@ref)             — read DEM raster
8. [`lift_to_3d`](@ref)           — sample boundary point elevations
9. [`generate_mesh_volume`](@ref) — flat 2D mesh → prism extrusion → 3 tets per prism
"""
module GeoGmsh

import ArchGDAL
using DataFrames
using GeoDataFrames
import GeoInterface as GI
import GeometryOps as GO
using Printf
import Proj

include("geometry.jl")
include("ingest.jl")
include("io.jl")
include("simplify.jl")
include("projection.jl")
include("terrain.jl")
include("adaptivity.jl")
include("tiles.jl")
include("verbose.jl")
include("gmsh.jl")
include("pipeline.jl")

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

export Geometry2D, Geometry3D, Contour, npoints, nedges

export ingest

export read_geodata, list_components, read_shapefile

export MinEdgeLength, AngleFilter, ComposedAlg

export rescale, filter_components

export AdaptivityAlgorithm, nominal_size, apply_adaptivity!
export SlopeAdaptivity, BoundaryLayerAdaptivity, ComposedAdaptivity

export DEMRaster, read_dem, sample_elevation, lift_to_3d

export download_cop30_tiles, prepare_dem

export write_geo
export generate_mesh
export generate_mesh_volume

export geoms_to_geo
export geoms_to_msh
export geoms_to_geo_3d
export geoms_to_msh_3d
export geoms_to_msh_3d_volume
export shapefile_to_geo
export shapefile_to_msh

end
