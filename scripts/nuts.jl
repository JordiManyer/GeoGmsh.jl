using ShapefileToGmsh

input  = joinpath(@__DIR__, "..", "meshes", "NUTS", "NUTS_RG_01M_2024_3035.shp")
output = joinpath(@__DIR__, "..", "data", "nuts")

comps = list_components(input);
sort!(filter(row -> row.NUTS_ID ∈ ("ES","ES51"), comps), by = row -> row.area, rev=true)[1:5]

# CRS: ETRS89 LAEA (EPSG:3035), already in metres — skip reprojection.
# Coarsen to 500 km minimum edge length before meshing.
# Rescale into a 100×100 bounding box so mesh_size is in those units.
shapefile_to_msh(
  input, output;
  select            = row -> row.NUTS_ID ∈ ("ES", "ES51") && row.idx ∈ (1756, 1490) && row.ring == 1,
  name_fn           = row -> string(row.NUTS_ID),
  proj_method       = nothing,
  edge_length_range = (10_000.0, Inf),
  bbox_size         = 100.0,
  mesh_size         = 2.0,
  split_components  = true,
)

println("Written: $output/")
