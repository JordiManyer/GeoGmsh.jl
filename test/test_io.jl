module IOTests

using ShapefileToGmsh
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")

function run()
  # read_shapefile: basic structure.
  geoms, crs = read_shapefile(AUS_SHP)
  @test !isempty(geoms)
  @test crs isa String   # raw WKT from .prj
  for g in geoms
    @test g.exterior.closed
    @test npoints(g.exterior) >= 3
    for h in g.holes
      @test h.closed
      @test npoints(h) >= 3
    end
  end

  # list_components: returns a NamedTuple vector with :idx.
  meta = list_components(AUS_SHP)
  @test !isempty(meta)
  @test haskey(meta[1], :idx)
  @test meta[1].idx == 1

  # read_shapefile with select by index.
  geoms_sel, _ = read_shapefile(AUS_SHP; select = [1])
  @test length(geoms_sel) == length(geoms)   # Australia has 1 record; same result

  # read_shapefile with predicate (include all — predicate always true).
  geoms_all, _ = read_shapefile(AUS_SHP; select = _ -> true)
  @test length(geoms_all) == length(geoms)

  # read_shapefile with predicate (exclude all).
  geoms_none, _ = read_shapefile(AUS_SHP; select = _ -> false)
  @test isempty(geoms_none)
end

end # module
