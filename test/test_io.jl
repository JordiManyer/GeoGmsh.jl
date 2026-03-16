module IOTests

using ShapefileToGmsh
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")

function run()
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
end

end # module
