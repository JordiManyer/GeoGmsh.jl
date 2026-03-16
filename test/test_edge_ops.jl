module EdgeOpsTests

using ShapefileToGmsh
using Test

# Dense circle polygon: N equispaced points on a circle of radius R.
# Adjacent points are ≈ 2πR/N apart.
function _dense_circle(n::Int, R::Float64 = 50_000.0)
  pts = NTuple{2,Float64}[
    (R * cos(2π * i / n), R * sin(2π * i / n)) for i in 0:n-1
  ]
  ShapeGeometry(Contour(pts, true), Contour[])
end

function run()
  # Two synthetic components: 400-point circle (≈785 m spacing) and
  # 200-point circle (≈628 m spacing).
  geoms = [_dense_circle(400, 50_000.0), _dense_circle(200, 20_000.0)]
  _test_coarsening(geoms)
  _test_refinement(geoms)
  _test_filter()
end

function _test_coarsening(geoms)
  total_before = sum(npoints(g.exterior) for g in geoms)

  coarsened_s = coarsen_edges(geoms, 5_000.0; strategy = :single)
  coarsened_i = coarsen_edges(geoms, 5_000.0; strategy = :iterative)

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

function _test_refinement(geoms)
  total_before = sum(npoints(g.exterior) for g in geoms)
  refined      = refine_edges(geoms, 1_000.0)
  total_after  = sum(npoints(g.exterior) for g in refined)
  @test total_after >= total_before
end

function _test_filter()
  # A 3-point triangle should be dropped by filter_components (min_points=4).
  tri = ShapeGeometry(
    Contour(NTuple{2,Float64}[(0.0,0.0),(1.0,0.0),(0.5,1.0)], true),
    Contour[],
  )
  square = _dense_circle(10, 1.0)   # 10-point ring: kept
  filtered = filter_components([tri, square])
  @test length(filtered) == 1
  @test npoints(filtered[1].exterior) == 10

  # min_points=3 should keep the triangle.
  @test length(filter_components([tri, square]; min_points = 3)) == 2

  # Degenerate hole removed, exterior kept.
  tri_hole = ShapeGeometry(
    square.exterior,
    [Contour(NTuple{2,Float64}[(0.0,0.0),(0.1,0.0),(0.05,0.1)], true)],
  )
  filtered2 = filter_components([tri_hole])
  @test length(filtered2) == 1
  @test isempty(filtered2[1].holes)
end

end # module
