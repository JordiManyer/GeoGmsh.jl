module ProjectionTests

using ShapefileToGmsh
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")

function run()
  geoms, crs = read_shapefile(AUS_SHP)

  # Symbol shorthand and struct types must both work.
  for method in (:equirectangular, :mercator, :web_mercator,
                  Equirectangular(), Mercator(), WebMercator(),
                  Equirectangular(lat_ref = -25.0))
    proj = project_to_meters(geoms, crs; method)
    @test length(proj) == length(geoms)
    # Australia spans roughly lon 113–154, lat -44 to -10 (degrees).
    # After projection the x-coordinates should be in the ~10⁷ m range.
    x_vals = [pt[1] for g in proj for pt in g.exterior.points]
    @test maximum(abs, x_vals) > 1e6
  end

  # Identity when CRS is already metres.
  geoms_m = project_to_meters(geoms, :meters)
  @test geoms_m === geoms
end

end # module
