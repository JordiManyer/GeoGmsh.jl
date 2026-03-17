module MinEdgeLengthTests

using GeoGmsh
import GeometryOps as GO
import GeoInterface as GI
using Test

# Build a polygon from an explicit list of (x,y) pairs.
function _poly(pts)
  closed = vcat(pts, [pts[1]])
  GI.Polygon([GI.LinearRing(closed)])
end

# Count points in the exterior ring (excluding the repeated closing point).
function _npts(poly)
  ring = GI.getexterior(poly)
  n = GI.npoint(ring)
  # GI rings include the closing point
  pts = [(GI.x(GI.getpoint(ring, i)), GI.y(GI.getpoint(ring, i))) for i in 1:n]
  first(pts) == last(pts) ? n - 1 : n
end

function run()
  _test_basic_reduction()
  _test_strict_min_edge_guarantee()
  _test_endpoints_preserved()
  _test_trivial_inputs()
  _test_only_tol_allowed()
end

function _test_basic_reduction()
  # 8 equispaced points on a unit circle → spacing ≈ 0.765.
  # With tol=0.9 all interior points should be removed except the endpoints.
  n = 8
  pts = [(cos(2π*i/n), sin(2π*i/n)) for i in 0:n-1]
  poly = _poly(pts)
  simplified = GO.simplify(MinEdgeLength(tol = 0.9), poly)
  @test _npts(simplified) < _npts(poly)
  @test _npts(simplified) >= 3
end

function _test_strict_min_edge_guarantee()
  # Place points at x = 0, 1, 1.4, 2.4, 2.5, 3.5 (y=0 line, closed as triangle).
  # After simplification with tol=1.1, no consecutive kept-pair should be < 1.1.
  pts = [(0.0,0.0),(1.0,0.0),(1.4,0.0),(2.4,0.0),(2.5,0.0),(3.5,1.0)]
  poly = _poly(pts)
  simplified = GO.simplify(MinEdgeLength(tol = 1.1), poly)
  ring = GI.getexterior(simplified)
  n = GI.npoint(ring)
  for i in 1:(n-1)
    p1 = GI.getpoint(ring, i)
    p2 = GI.getpoint(ring, i+1)
    d  = hypot(GI.x(p2)-GI.x(p1), GI.y(p2)-GI.y(p1))
    # Every edge except possibly the closing wrap-around must be >= tol.
    @test d >= 1.1 - 1e-10
  end
end

function _test_endpoints_preserved()
  # The first and last points of the ring must always be retained.
  n = 20
  pts = [(cos(2π*i/n), sin(2π*i/n)) for i in 0:n-1]
  poly = _poly(pts)
  simplified = GO.simplify(MinEdgeLength(tol = 0.5), poly)
  ring_in  = GI.getexterior(poly)
  ring_out = GI.getexterior(simplified)
  p_first_in  = GI.getpoint(ring_in,  1)
  p_first_out = GI.getpoint(ring_out, 1)
  @test GI.x(p_first_in) ≈ GI.x(p_first_out)
  @test GI.y(p_first_in) ≈ GI.y(p_first_out)
end

function _test_trivial_inputs()
  # Two-point degenerate input should pass through unchanged.
  pts = [(0.0,0.0),(1.0,0.0)]
  ring = GI.LinearRing(vcat(pts, [pts[1]]))
  poly = GI.Polygon([ring])
  result = GO.simplify(MinEdgeLength(tol = 0.1), poly)
  @test !isnothing(result)
end

function _test_only_tol_allowed()
  # MinEdgeLength requires exactly one of tol/number/ratio.
  @test_throws Exception MinEdgeLength()                     # none set
  @test_throws Exception MinEdgeLength(tol=1.0, number=3)    # two set
  @test MinEdgeLength(tol=1.0) isa MinEdgeLength
end

end # module
