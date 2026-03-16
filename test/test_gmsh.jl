module GmshTests

using ShapefileToGmsh
using Test

const FIXTURE_SHP = joinpath(@__DIR__, "..", "meshes", "toy", "fixture.shp")
const TMP_NAME    = joinpath(tempdir(), "test_shapefile_to_gmsh")

function run()
  geoms, _ = read_shapefile(FIXTURE_SHP)

  # Single-file output → TMP_NAME.geo
  write_geo(geoms, TMP_NAME; mesh_size = 1.0)
  @test isfile(TMP_NAME * ".geo")
  content = read(TMP_NAME * ".geo", String)
  @test occursin("Point(",         content)
  @test occursin("Line(",          content)
  @test occursin("Curve Loop(",    content)
  @test occursin("Plane Surface(", content)
  rm(TMP_NAME * ".geo")

  # Split-components output → TMP_NAME/ directory, one file per geometry.
  write_geo(geoms, TMP_NAME; mesh_size = 1.0, split_components = true)
  @test isdir(TMP_NAME)
  @test length(readdir(TMP_NAME)) == length(geoms)
  @test all(endswith(".geo"), readdir(TMP_NAME))
  rm(TMP_NAME; recursive = true)
end

end # module
