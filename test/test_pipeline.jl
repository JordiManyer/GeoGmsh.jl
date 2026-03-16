module PipelineTests

using ShapefileToGmsh
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")
const TMP_NAME = joinpath(tempdir(), "test_shapefile_to_gmsh_pipeline")

function run()
  # Target CRS string, edge coarsening, and bbox rescaling.
  shapefile_to_geo(
    AUS_SHP, TMP_NAME;
    proj_method       = "EPSG:3857",
    edge_length_range = (500_000.0, Inf),
    bbox_size         = 100.0,
    mesh_size         = 2.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # proj_method = nothing (skip reprojection, keep raw coordinates).
  shapefile_to_geo(
    AUS_SHP, TMP_NAME;
    proj_method = nothing,
    mesh_size   = 1.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")
end

end # module
