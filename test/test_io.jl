module IOTests

using ShapefileToGmsh
using Test

include("fixture.jl")
const FIXTURE_SHP = _create_test_shapefile(mktempdir())

function run()
  # read_shapefile: basic structure (3 records → 3 ShapeGeometry objects).
  geoms, crs = read_shapefile(FIXTURE_SHP)
  @test length(geoms) == 3
  @test crs isa String
  @test occursin("GEOGCS", crs)
  for g in geoms
    @test g.exterior.closed
    @test npoints(g.exterior) >= 3
  end

  # list_components: returns one NamedTuple per record with :idx and DBF columns.
  meta = list_components(FIXTURE_SHP)
  @test length(meta) == 3
  @test meta[1].idx  == 1
  @test string(meta[1].NAME) == "Alpha"
  @test string(meta[2].CODE) == "B2"
  @test Int(meta[3].VAL)     == 3

  # select by index vector.
  g12, _ = read_shapefile(FIXTURE_SHP; select = [1, 2])
  @test length(g12) == 2

  g3, _ = read_shapefile(FIXTURE_SHP; select = [3])
  @test length(g3) == 1

  # select by predicate.
  gb, _ = read_shapefile(FIXTURE_SHP; select = row -> string(row.CODE) == "B2")
  @test length(gb) == 1

  # select = nothing → all records.
  gall, _ = read_shapefile(FIXTURE_SHP; select = nothing)
  @test length(gall) == 3

  # select with predicate that excludes all.
  gnone, _ = read_shapefile(FIXTURE_SHP; select = _ -> false)
  @test isempty(gnone)
end

end # module
