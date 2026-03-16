# ShapefileToGmsh.jl

[![Build Status](https://github.com/JordiManyer/ShapefileToGmsh.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JordiManyer/ShapefileToGmsh.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JordiManyer.github.io/ShapefileToGmsh.jl/dev/)

A Julia package that converts ESRI Shapefiles into Gmsh geometry (`.geo`) and mesh (`.msh`) files, with full control over coordinate projection, edge resolution, and component selection.

<table>
  <tr>
    <td align="center"><b>Australia — geometry</b></td>
    <td align="center"><b>Australia — mesh</b></td>
  </tr>
  <tr>
    <td><img src="docs/assets/australia_geo.png" alt="Australia geometry" width="100%"/></td>
    <td><img src="docs/assets/australia_mesh.png" alt="Australia mesh" width="100%"/></td>
  </tr>
  <tr>
    <td align="center"><b>Spain</b></td>
    <td align="center"><b>Catalonia</b></td>
  </tr>
  <tr>
    <td><img src="docs/assets/spain.png" alt="Spain mesh" width="100%"/></td>
    <td><img src="docs/assets/catalonia.png" alt="Catalonia mesh" width="100%"/></td>
  </tr>
</table>

## Features

- **Read & filter** — load Shapefiles and inspect their attribute tables; filter records or individual rings by any DBF attribute, geometry size, or bounding box.
- **Reproject** — convert between coordinate systems via [Proj.jl](https://github.com/JuliaGeo/Proj.jl) (e.g. geographic degrees → Web Mercator metres).
- **Resample edges** — coarsen over-resolved coastlines or refine coarse boundaries to a target edge length.
- **Rescale** — normalise geometry into a dimensionless bounding box so `mesh_size` stays meaningful across datasets.
- **Output** — write a human-readable `.geo` script (open in the Gmsh GUI) or call the Gmsh API directly to produce a `.msh` file; linear or quadratic elements, triangle or quad recombination.
- **Split components** — one file per polygon ring, with user-defined filenames via `name_fn`.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/JordiManyer/ShapefileToGmsh.jl")
```

## Quick start

```julia
using ShapefileToGmsh

# 1. Inspect the shapefile — see what records and rings are available.
list_components("NUTS_RG_01M_2024_3035.shp")

# 2. Mesh Spain and Catalonia, one file each, named by NUTS ID.
shapefile_to_msh(
  "NUTS_RG_01M_2024_3035.shp",
  "output/nuts";
  select            = row -> row.NUTS_ID ∈ ("ES", "ES51") && row.ring == 1,
  name_fn           = row -> string(row.NUTS_ID),
  proj_method       = nothing,          # already in metres (EPSG:3035)
  edge_length_range = (10_000.0, Inf),  # coarsen to ≥ 10 km edges
  bbox_size         = 100.0,
  mesh_size         = 2.0,
  split_components  = true,
)
# → output/nuts/ES.msh, output/nuts/ES51.msh
```

See the [documentation](https://JordiManyer.github.io/ShapefileToGmsh.jl/dev/) for the full pipeline guide and API reference.

## Data sources

The example meshes above were generated from the following open datasets:

- **NUTS — Territorial units statistics** (Spain, Catalonia): Eurostat / GISCO, © European Union, [https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics](https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics)
- **ASGS Edition 3** (Australia): Australian Bureau of Statistics, [https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3](https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3)
