module PipelineTests

using ShapefileToGmsh
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")
const TMP_NAME = joinpath(tempdir(), "test_shapefile_to_gmsh_pipeline")

function run()
  # Default projection (equirectangular symbol), edge range, single file.
  shapefile_to_geo(
    AUS_SHP, TMP_NAME;
    edge_length_range = (50_000.0, Inf),
    mesh_size         = 50_000.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # Struct-based projection method.
  shapefile_to_geo(
    AUS_SHP, TMP_NAME;
    proj_method       = Equirectangular(lat_ref = -25.0),
    edge_length_range = (50_000.0, Inf),
    mesh_size         = 50_000.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")
end

end # module
