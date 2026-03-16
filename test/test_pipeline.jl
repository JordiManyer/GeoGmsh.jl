module PipelineTests

using ShapefileToGmsh
using Test

include("fixture.jl")
const FIXTURE_SHP = _create_test_shapefile(mktempdir())
const TMP_NAME    = joinpath(tempdir(), "test_shapefile_to_gmsh_pipeline")

function run()
  # Reproject + coarsen (edge_length_range much larger than the 2° squares,
  # so rings survive) + rescale → single .geo file.
  shapefile_to_geo(
    FIXTURE_SHP, TMP_NAME;
    proj_method       = "EPSG:3857",
    edge_length_range = (50_000.0, Inf),
    bbox_size         = 100.0,
    mesh_size         = 5.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # Skip reprojection.
  shapefile_to_geo(FIXTURE_SHP, TMP_NAME; proj_method = nothing, mesh_size = 1.0)
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # select kwarg filters records before processing.
  shapefile_to_geo(
    FIXTURE_SHP, TMP_NAME;
    proj_method = nothing,
    select      = row -> string(row.CODE) == "A1",
    mesh_size   = 1.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")
end

end # module
