using ShapefileToGmsh

input  = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")
output = joinpath(@__DIR__, "..", "data", "australia")

comps = list_components(input);
sort!(comps, by = row -> row.area, rev=true)[1:5]

shapefile_to_geo(
  input, output;
  select            = row -> (row.AUS_CODE21 == "AUS" && row.ring == 1),
  name_fn           = row -> string(row.AUS_CODE21),
  proj_method       = "EPSG:3857",
  edge_length_range = (100_000.0, Inf),
  bbox_size         = 100.0,
  verbose           = true,
  split_components  = true,
)

# CRS: GDA2020 geographic (degrees) → reproject to Web Mercator (metres).
# Coarsen to 500 km minimum edge length before meshing.
# Rescale into a 100×100 bounding box so mesh_size is in those units.
shapefile_to_msh(
  input, output;
  select            = row -> (row.AUS_CODE21 == "AUS" && row.ring == 1),
  name_fn           = row -> string(row.AUS_CODE21),
  proj_method       = "EPSG:3857",
  edge_length_range = (100_000.0, Inf),
  bbox_size         = 100.0,
  mesh_size         = 2.0,
  verbose           = true,
  split_components  = true,
)
