using ShapefileToGmsh

input  = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")
output = joinpath(@__DIR__, "..", "data", "australia")

# CRS: GDA2020 geographic (degrees) → reproject to Web Mercator (metres).
# Coarsen to 500 km minimum edge length before meshing.
# Rescale into a 100×100 bounding box so mesh_size is in those units.
shapefile_to_msh(
  input, output;
  proj_method       = "EPSG:3857",
  edge_length_range = (500_000.0, Inf),
  bbox_size         = 100.0,
  mesh_size         = 2.0,
  split_components  = true,
)

println("Written: $output/")
