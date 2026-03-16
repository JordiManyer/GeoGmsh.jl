module PipelineTests

using GeoGmsh
import GeoInterface as GI
using Test

const FIXTURE_SHP = joinpath(@__DIR__, "..", "meshes", "toy", "fixture.shp")
const TMP_NAME    = joinpath(tempdir(), "test_geogmsh_pipeline")

function run()
  # geoms_to_geo from a DataFrame: reproject + simplify + rescale
  df = read_geodata(FIXTURE_SHP)
  geoms_to_geo(df, TMP_NAME;
    target_crs   = "EPSG:3857",
    simplify_tol = 50_000.0,
    bbox_size    = 100.0,
    mesh_size    = 5.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # Skip reprojection
  geoms_to_geo(df, TMP_NAME; target_crs = nothing, mesh_size = 1.0)
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # select kwarg filters rows before processing
  geoms_to_geo(df, TMP_NAME;
    target_crs = nothing,
    select     = row -> string(row.NAME) == "Alpha",
    mesh_size  = 1.0,
  )
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # backward-compat wrapper
  shapefile_to_geo(FIXTURE_SHP, TMP_NAME;
    target_crs = nothing, mesh_size = 1.0, verbose = false)
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")

  # geoms_to_geo from raw GeoInterface geometry (no CRS)
  poly = GI.Polygon([GI.LinearRing([(0.,0.),(1.,0.),(1.,1.),(0.,1.),(0.,0.)])])
  geoms_to_geo(poly, TMP_NAME; target_crs = nothing, mesh_size = 0.1)
  @test isfile(TMP_NAME * ".geo")
  rm(TMP_NAME * ".geo")
end

end # module
