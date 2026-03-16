module GmshTests

using ShapefileToGmsh
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")
const TMP_NAME = joinpath(tempdir(), "test_shapefile_to_gmsh")

function run()
  geoms, crs = read_shapefile(AUS_SHP)
  proj      = project_to_meters(geoms, crs)
  coarsened = coarsen_edges(proj, 50_000.0)

  # Single-file output → TMP_NAME.geo
  write_geo(coarsened, TMP_NAME; mesh_size = 50_000.0)
  @test isfile(TMP_NAME * ".geo")
  content = read(TMP_NAME * ".geo", String)
  @test occursin("Point(",         content)
  @test occursin("Line(",          content)
  @test occursin("Curve Loop(",    content)
  @test occursin("Plane Surface(", content)
  rm(TMP_NAME * ".geo")

  # Split-components output → TMP_NAME/ directory
  write_geo(coarsened, TMP_NAME; mesh_size = 50_000.0, split_components = true)
  @test isdir(TMP_NAME)
  component_files = readdir(TMP_NAME)
  @test length(component_files) == length(coarsened)
  @test all(endswith(".geo"), component_files)
  rm(TMP_NAME; recursive = true)
end

end # module
