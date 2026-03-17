# GeoGmsh.jl

[![Build Status](https://github.com/JordiManyer/GeoGmsh.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JordiManyer/GeoGmsh.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JordiManyer.github.io/GeoGmsh.jl/dev/)

A Julia package that converts geospatial data into Gmsh geometry (`.geo`) and
mesh (`.msh`) files. Accepts any GeoInterface-compatible source: Shapefiles,
GeoJSON, GeoPackage, GeoParquet, NaturalEarth data, or raw geometries.

**Flat 2D meshes** from any geospatial boundary:

<table>
  <tr>
    <td align="center"><b>Iberian Peninsula</b></td>
    <td align="center"><b>Catalonia</b></td>
    <td align="center"><b>Australia — geometry</b></td>
    <td align="center"><b>Australia — mesh</b></td>
  </tr>
  <tr>
    <td><img src="docs/src/assets/iberia.png" alt="Iberia" width="100%"/></td>
    <td><img src="docs/src/assets/catalonia.png" alt="Catalonia" width="100%"/></td>
    <td><img src="docs/src/assets/australia_geo.png" alt="Australia geometry" width="100%"/></td>
    <td><img src="docs/src/assets/australia_mesh.png" alt="Australia mesh" width="100%"/></td>
  </tr>
</table>

**3D terrain manifolds** by sampling a Digital Elevation Model at every mesh node:

<table>
  <tr>
    <td align="center"><b>Mont Blanc</b></td>
  </tr>
  <tr>
    <td><img src="docs/src/assets/montblanc.png" alt="Mont Blanc terrain mesh" width="60%"/></td>
  </tr>
</table>

## Features

- **Universal reader** — load any geospatial format (Shapefile, GeoJSON, GeoPackage, GeoParquet, …) through a single `read_geodata` call backed by GDAL. Read directly from ZIP archives via `/vsizip/`.
- **Reproject** — convert between coordinate systems via Proj.jl (e.g. geographic degrees → UTM metres).
- **Simplify** — `MinEdgeLength` removes short edges; `AngleFilter` suppresses zig-zag spikes; compose algorithms with `∘`.
- **Segmentize** — subdivide long edges to control maximum element size.
- **Rescale** — normalise geometry into a dimensionless bounding box so `mesh_size` is consistent.
- **3D terrain** — lift a 2D mesh to terrain elevation by sampling a DEM raster at every node.
- **Output** — `.geo` script or `.msh` file via the Gmsh API; linear/quadratic elements, triangle or quad recombination.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/JordiManyer/GeoGmsh.jl")
```

## Quick start

```julia
using GeoGmsh, NaturalEarth

# No files needed — NaturalEarth.jl provides built-in country boundaries
countries = naturalearth("admin_0_countries", 110)

geoms_to_msh(countries, "france";
  select       = row -> get(row, :NAME, "") == "France" && row.ring == 1,
  target_crs   = "EPSG:3857",
  simplify_alg = MinEdgeLength(tol = 5_000.0),
  bbox_size    = 100.0,
  mesh_size    = 2.0,
)
# → france.msh
```

```julia
using GeoGmsh

# 3D terrain mesh from a bounding box + DEM
geoms_to_msh_3d("bbox.geojson", "dem_utm.tif", "terrain";
  target_crs   = "EPSG:32632",
  simplify_alg = MinEdgeLength(tol = 500.0),
  mesh_size    = 500.0,
)
# → terrain.msh  (node z-coordinates sampled from DEM)
```

See the [documentation](https://JordiManyer.github.io/GeoGmsh.jl/dev/) for the
full pipeline guide, API reference, and worked examples.

## Data sources

- **Natural Earth**: [naturalearthdata.com](https://www.naturalearthdata.com/) — free vector map data.
- **NUTS**: Eurostat / GISCO, © European Union.
- **ASGS Edition 3**: Australian Bureau of Statistics.
- **Copernicus GLO-30 DEM**: © DLR / Airbus.
