module ProjectionTests

using ShapefileToGmsh
import Proj
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")

function run()
  geoms, source_crs = read_shapefile(AUS_SHP)
  @test source_crs isa String   # raw WKT from .prj

  # Project using a target CRS string.
  proj1 = project_to_meters(geoms, source_crs; target = "EPSG:3857")
  @test length(proj1) == length(geoms)
  # Australia spans roughly lon 113–154. After Web Mercator the x-coordinates
  # should be in the ~10⁷ m range.
  x_vals = [pt[1] for g in proj1 for pt in g.exterior.points]
  @test maximum(abs, x_vals) > 1e6

  # Project using a pre-built Proj.Transformation.
  trans = Proj.Transformation(source_crs, "EPSG:3857"; always_xy = true)
  proj2 = project_to_meters(geoms, source_crs; target = trans)
  @test length(proj2) == length(geoms)

  # rescale into a 100×100 bounding box.
  scaled = rescale(proj1, 100.0)
  xs = [pt[1] for g in scaled for pt in g.exterior.points]
  ys = [pt[2] for g in scaled for pt in g.exterior.points]
  @test minimum(xs) ≈ 0.0  atol = 1e-8
  @test minimum(ys) ≈ 0.0  atol = 1e-8
  @test max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)) ≈ 100.0  atol = 1e-6
end

end # module
