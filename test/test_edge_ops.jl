module EdgeOpsTests

using ShapefileToGmsh
using Test

const AUS_SHP = joinpath(@__DIR__, "..", "meshes", "australia", "AUS_2021_AUST_GDA2020.shp")

function run()
  geoms, crs = read_shapefile(AUS_SHP)
  proj = project_to_meters(geoms, crs)

  _test_coarsening(proj)
  _test_refinement(proj)
end

function _test_coarsening(proj)
  total_before = sum(npoints(g.exterior) for g in proj)

  coarsened_s = coarsen_edges(proj, 50_000.0; strategy=:single)
  coarsened_i = coarsen_edges(proj, 50_000.0; strategy=:iterative)

  total_single    = sum(npoints(g.exterior) for g in coarsened_s)
  total_iterative = sum(npoints(g.exterior) for g in coarsened_i)

  @test total_single    < total_before
  @test total_iterative < total_before
  # Iterative must be at least as aggressive as single-pass.
  @test total_iterative <= total_single

  # No ring should degenerate below 3 points.
  for g in coarsened_i
    @test npoints(g.exterior) >= 3
  end
end

function _test_refinement(proj)
  total_before = sum(npoints(g.exterior) for g in proj)
  refined      = refine_edges(proj, 10_000.0)   # 10 km max edge
  total_after  = sum(npoints(g.exterior) for g in refined)
  @test total_after >= total_before
end

end # module
