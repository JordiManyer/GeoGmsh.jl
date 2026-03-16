module IOTests

using GeoGmsh
using DataFrames
import GeoInterface as GI
using Test

const FIXTURE_SHP = joinpath(@__DIR__, "..", "meshes", "toy", "fixture.shp")

function run()
  # read_geodata: basic structure
  df = read_geodata(FIXTURE_SHP)
  @test nrow(df) == 3
  @test :geometry in propertynames(df)
  @test :NAME in propertynames(df)

  # CRS should be present (WGS84 from the .prj sidecar)
  crs = GI.crs(df)
  @test !isnothing(crs)
  @test occursin("4326", string(crs.val))

  # All geometries should be valid polygons
  for row in eachrow(df)
    @test !isnothing(row.geometry)
    @test GI.geomtrait(row.geometry) isa GI.PolygonTrait
  end

  # select by predicate
  df2 = read_geodata(FIXTURE_SHP; select = row -> string(row.NAME) == "Beta")
  @test nrow(df2) == 1
  @test string(df2.NAME[1]) == "Beta"

  # select = nothing → all records
  df_all = read_geodata(FIXTURE_SHP; select = nothing)
  @test nrow(df_all) == 3

  # select that excludes all
  df_none = read_geodata(FIXTURE_SHP; select = _ -> false)
  @test nrow(df_none) == 0

  # list_components: returns a DataFrame with geometry stats
  df_listed = list_components(FIXTURE_SHP)
  @test nrow(df_listed) == 3
  @test :n_pts in propertynames(df_listed)
  @test :area  in propertynames(df_listed)
  @test :xmin  in propertynames(df_listed)
  @test all(df_listed.n_pts .> 0)

  # read_shapefile is a backward-compatible alias
  df3 = read_shapefile(FIXTURE_SHP)
  @test nrow(df3) == nrow(df)
end

end # module
