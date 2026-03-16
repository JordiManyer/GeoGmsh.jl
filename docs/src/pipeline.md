# [Pipeline guide](@id pipeline)

This page walks through the full conversion pipeline step by step.
All steps are also available individually so you can mix and match.

## 1. Inspect the attribute table

Official shapefiles often contain many components — countries, statistical
regions, islands — identified by attribute columns in the `.dbf` sidecar.
Call [`list_components`](@ref) before loading geometry to understand what is
available and find the field values to filter on:

```julia
meta = list_components("NUTS_RG_01M_2024_3035.shp")
# prints a formatted table with idx, NUTS_ID, LEVL_CODE, CNTR_CODE, NAME_LATN, …

# Access programmatically:
germany = filter(r -> r.CNTR_CODE == "DE", meta)
```

The returned `Vector{NamedTuple}` has an `:idx` field (1-based row number) and
one field per DBF column.

## 2. Read the shapefile

```julia
geoms, source_crs = read_shapefile("file.shp")
```

`geoms` is a `Vector{ShapeGeometry}`.  Each `ShapeGeometry` holds one exterior
polygon ring and zero or more hole rings.  MultiPolygon records are
automatically flattened: each outer ring becomes its own `ShapeGeometry`.

`source_crs` is the raw WKT string from the `.prj` sidecar (or `nothing` if no
`.prj` is present).  It is passed directly to PROJ for reprojection.

### Selecting components

Pass `select` to load only a subset of records:

```julia
# By row index (as shown by list_components)
geoms, crs = read_shapefile("file.shp"; select = [3, 7, 12])

# By predicate on DBF attributes
geoms, crs = read_shapefile("NUTS.shp";
  select = row -> row.CNTR_CODE == "DE" && row.LEVL_CODE == 0)
```

!!! tip
    Use `list_components` first to discover field names and values.

## 3. Reproject coordinates

Geographic shapefiles store coordinates in decimal degrees (lon, lat).
Gmsh works in a flat Cartesian plane, so reprojection to metres is required
before edge-length operations or meshing.

```julia
geoms = project_to_meters(geoms, source_crs; target = "EPSG:3857")
```

`target` may be any CRS string understood by PROJ (EPSG code, WKT, PROJ
string) or a pre-built `Proj.Transformation`:

```julia
import Proj
trans = Proj.Transformation(source_crs, "EPSG:3035"; always_xy = true)
geoms = project_to_meters(geoms, source_crs; target = trans)
```

Pass `proj_method = nothing` to the pipeline functions (or skip this step
entirely) when the shapefile is already in metres — for example NUTS data in
EPSG:3035.

## 4. Adjust edge resolution

Large shapefiles contain far more points than a coarse mesh needs, and
generating the mesh from millions of points is slow.

### Coarsening

Remove intermediate points to enforce a minimum edge length:

```julia
geoms = coarsen_edges(geoms, 500_000.0)   # ≥ 500 km edges
```

Two strategies are available via the `strategy` keyword:

| Strategy | Behaviour |
|----------|-----------|
| `:iterative` (default) | Repeat single-pass sweeps until stable. Guarantees no edge shorter than the threshold remains. |
| `:single` | One left-to-right sweep. Fast; may leave residual short edges. |

### Refinement

Subdivide edges that are too long:

```julia
geoms = refine_edges(geoms, 100_000.0)   # ≤ 100 km edges
```

Edges longer than `max_edge_length` are split into equispaced sub-edges.

### Filtering degenerate components

Aggressive coarsening can reduce small islands to triangular rings (3 points)
that cause Gmsh to fail.  [`filter_components`](@ref) removes them:

```julia
geoms = filter_components(geoms)          # drops rings with < 4 points
geoms = filter_components(geoms; min_points = 5)   # stricter threshold
```

The pipeline functions call this automatically after coarsening.

### Combined edge-length range

In the pipeline functions, pass `edge_length_range = (min, max)` to run
coarsening and refinement together:

```julia
shapefile_to_geo(...; edge_length_range = (50_000.0, 500_000.0))
```

## 5. Rescale geometry

When mesh parameters are easier to express in normalised units (e.g. for an
FEM solver that expects an O(1) domain), rescale the geometry so its largest
bounding-box dimension equals `L`:

```julia
geoms = rescale(geoms, 100.0)   # fits into a 100 × 100 box, origin at (0,0)
```

After rescaling, set `mesh_size` in the same normalised units (e.g. `2.0`
gives roughly 50 elements across the domain).

## 6. Write output

### `.geo` script

Produces a human-readable Gmsh geometry script that you can open in the GUI,
edit manually, or mesh with `gmsh file.geo -2`:

```julia
write_geo(geoms, "output/germany"; mesh_size = 2.0)
# → output/germany.geo
```

### `.msh` file (via Gmsh API)

Generates the mesh immediately without any intermediate file:

```julia
generate_mesh(geoms, "output/germany";
  mesh_size = 2.0,
  order     = 2,        # quadratic elements
  recombine = true,     # triangles → quads
)
# → output/germany.msh
```

### Split components

Both output functions accept `split_components = true` to write one file per
`ShapeGeometry` into a directory:

```julia
generate_mesh(geoms, "output/australia"; mesh_size = 2.0, split_components = true)
# → output/australia/0001.msh, 0002.msh, …
```

## One-call pipeline

For most use cases `shapefile_to_geo` or `shapefile_to_msh` compose the
entire pipeline:

```julia
shapefile_to_msh(
  "NUTS_RG_01M_2024_3035.shp",
  "output/germany";
  select            = row -> row.CNTR_CODE == "DE" && row.LEVL_CODE == 0,
  proj_method       = nothing,
  edge_length_range = (50_000.0, Inf),
  bbox_size         = 100.0,
  mesh_size         = 2.0,
  split_components  = true,
)
```

## Shapefile conventions

### Winding order

ESRI Shapefiles use the **opposite** winding convention to the OGC standard:

- Clockwise (CW) rings → exterior polygons
- Counter-clockwise (CCW) rings → holes

`read_shapefile` detects this automatically using the sign of the shoelace
area and normalises ring orientation for Gmsh:

- Stored exteriors → CCW (positive signed area)
- Stored holes → CW (negative signed area)

### MultiPolygon flattening

A single Shapefile record (one row in the attribute table) may contain
multiple disjoint outer rings — an archipelago stored as one feature, for
example.  Each outer ring (plus its holes) is returned as a separate
`ShapeGeometry`.  Consequently, `length(geoms)` may exceed the number of
DBF records.
