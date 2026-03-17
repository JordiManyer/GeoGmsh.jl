# [API reference](@id api)

## Pipeline

```@docs
geoms_to_geo
geoms_to_msh
shapefile_to_geo
shapefile_to_msh
```

## I/O

```@docs
read_geodata
list_components
read_shapefile
```

## Simplification

```@docs
MinEdgeLength
AngleFilter
ComposedAlg
```

## Mesh adaptivity

```@docs
AdaptivityAlgorithm
nominal_size
apply_adaptivity!
SlopeAdaptivity
BoundaryLayerAdaptivity
ComposedAdaptivity
```

## Rescaling and filtering

```@docs
rescale
filter_components
```

## Ingest

```@docs
ingest
```

## Gmsh output

```@docs
write_geo
generate_mesh
```

## DEM tiles

```@docs
download_cop30_tiles
download_opentopo_dem
prepare_dem
```

## 3D terrain

```@docs
DEMRaster
read_dem
sample_elevation
lift_to_3d
geoms_to_geo_3d
geoms_to_msh_3d
```

## 3D volume

```@docs
geoms_to_msh_3d_volume
generate_mesh_volume
```

## Geometry types

```@docs
Geometry2D
Geometry3D
Contour
npoints
nedges
```
