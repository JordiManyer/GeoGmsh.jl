module EdgeOpsTests

using GeoGmsh
import GeometryOps as GO
import GeoInterface as GI
using Test

# Dense circle as a GI.Polygon (for GO.simplify / GO.segmentize).
function _dense_circle_poly(n::Int, R::Float64 = 50_000.0)
  pts = [(R * cos(2π * i / n), R * sin(2π * i / n)) for i in 0:n-1]
  push!(pts, pts[1])   # close the ring
  GI.Polygon([GI.LinearRing(pts)])
end

# Dense circle as Geometry2D (for filter_components tests).
function _dense_circle_shape(n::Int, R::Float64 = 1.0)
  pts = NTuple{2,Float64}[(R * cos(2π * i / n), R * sin(2π * i / n)) for i in 0:n-1]
  Geometry2D(Contour(pts, true), Contour[])
end

function run()
  _test_simplify()
  _test_segmentize()
  _test_filter()
end

function _test_simplify()
  # 400-point circle with R=50 km → adjacent points ≈ 785 m apart.
  poly     = _dense_circle_poly(400, 50_000.0)
  n_before = GI.npoint(GI.getexterior(poly))

  simplified = GO.simplify(MinEdgeLength(tol = 5_000.0), poly)
  n_after    = GI.npoint(GI.getexterior(simplified))

  # Simplification must reduce the point count.
  @test n_after < n_before
  # Result must still be a valid ring (≥ 3 points).
  @test n_after >= 3
end

function _test_segmentize()
  # Coarse 4-point square-ish polygon → segmentize should add points.
  poly     = _dense_circle_poly(4, 50_000.0)
  n_before = GI.npoint(GI.getexterior(poly))

  segmented = GO.segmentize(poly; max_distance = 5_000.0)
  n_after   = GI.npoint(GI.getexterior(segmented))

  @test n_after > n_before
end

function _test_filter()
  # A 3-point triangle should be dropped by filter_components (min_points=4).
  tri = Geometry2D(
    Contour(NTuple{2,Float64}[(0.0,0.0),(1.0,0.0),(0.5,1.0)], true),
    Contour[],
  )
  square   = _dense_circle_shape(10, 1.0)
  filtered = filter_components([tri, square])
  @test length(filtered) == 1
  @test npoints(filtered[1].exterior) == 10

  # min_points=3 should keep the triangle.
  @test length(filter_components([tri, square]; min_points = 3)) == 2

  # Degenerate hole removed, exterior kept.
  tri_hole = Geometry2D(
    square.exterior,
    [Contour(NTuple{2,Float64}[(0.0,0.0),(0.1,0.0),(0.05,0.1)], true)],
  )
  filtered2 = filter_components([tri_hole])
  @test length(filtered2) == 1
  @test isempty(filtered2[1].holes)
end

end # module
