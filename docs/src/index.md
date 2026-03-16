# ShapefileToGmsh.jl

```@docs
ShapefileToGmsh
```

A Julia package that converts ESRI Shapefiles into Gmsh geometry (`.geo`) and
mesh (`.msh`) files, with full control over coordinate projection, edge
resolution, and component selection.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/you/ShapefileToGmsh.jl")   # or Pkg.develop
```

## Quick start

### Inspect a shapefile

Before meshing, inspect the attribute table to understand what components the
file contains and find the field values you want to filter on:

```julia
using ShapefileToGmsh

list_components("path/to/NUTS_RG_01M_2024_3035.shp")
```

```
idx   NUTS_ID   LEVL_CODE  CNTR_CODE  NAME_LATN
─────────────────────────────────────────────────
1     AL011     3          AL         Dibër
2     AL012     3          AL         Durrës
3     DE         0          DE         Deutschland
...
```

### Generate a `.geo` file

```julia
shapefile_to_geo(
  "NUTS_RG_01M_2024_3035.shp",
  "output/germany";
  select            = row -> row.CNTR_CODE == "DE" && row.LEVL_CODE == 0,
  proj_method       = nothing,          # already in metres (EPSG:3035)
  edge_length_range = (50_000.0, Inf),  # coarsen to ≥ 50 km edges
  bbox_size         = 100.0,            # normalise to 100 × 100
  mesh_size         = 2.0,
)
# → output/germany.geo
```

### Generate a `.msh` file directly

Replace `shapefile_to_geo` with `shapefile_to_msh` to invoke the Gmsh API and
produce a mesh in one step:

```julia
shapefile_to_msh(
  "AUS_2021_AUST_GDA2020.shp",
  "output/australia";
  proj_method       = "EPSG:3857",
  edge_length_range = (500_000.0, Inf),
  bbox_size         = 100.0,
  mesh_size         = 2.0,
  split_components  = true,   # one .msh per island
)
# → output/australia/0001.msh, 0002.msh, …
```

## Overview

The pipeline runs these steps in order:

| Step | Function | Purpose |
|------|----------|---------|
| Read | [`read_shapefile`](@ref) | Parse `.shp` geometry, filter records |
| Reproject | [`project_to_meters`](@ref) | Convert lon/lat → metres via PROJ |
| Coarsen | [`coarsen_edges`](@ref) | Remove points on short edges |
| Refine | [`refine_edges`](@ref) | Subdivide long edges |
| Filter | [`filter_components`](@ref) | Drop degenerate rings after coarsening |
| Rescale | [`rescale`](@ref) | Normalise into an L × L box |
| Output | [`write_geo`](@ref) / [`generate_mesh`](@ref) | Write `.geo` or `.msh` |

Each step is also available as a standalone function for more control.
See the [Pipeline guide](@ref pipeline) and [API reference](@ref api) for details.
