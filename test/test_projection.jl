module ProjectionTests

using ShapefileToGmsh
import Proj
using Test

# Small synthetic geometry at known geographic coordinates (degrees).
# A square near (lon=10°, lat=20°) with 1° side ≈ 111 km.
const _PTS = NTuple{2,Float64}[
  (10.0, 20.0), (11.0, 20.0), (11.0, 21.0), (10.0, 21.0),
]
const _GEOMS = [ShapeGeometry(Contour(_PTS, true), Contour[])]

function run()
  # Project to EPSG:3857 using a source CRS string.
  proj = project_to_meters(_GEOMS, "EPSG:4326"; target = "EPSG:3857")
  @test length(proj) == 1
  # At lon≈10°, x ≈ R·deg2rad(10) ≈ 1.11e6 m.
  x_vals = [pt[1] for pt in proj[1].exterior.points]
  @test all(x -> 1.0e6 < x < 2.0e6, x_vals)

  # Project using a pre-built Proj.Transformation.
  trans = Proj.Transformation("EPSG:4326", "EPSG:3857"; always_xy = true)
  proj2 = project_to_meters(_GEOMS, "EPSG:4326"; target = trans)
  @test length(proj2) == 1
  x2 = [pt[1] for pt in proj2[1].exterior.points]
  @test x2 ≈ x_vals   atol = 1.0

  # nothing source_crs warns and falls back to EPSG:4326.
  proj3 = (@test_logs (:warn,) project_to_meters(_GEOMS, nothing; target = "EPSG:3857"))
  @test length(proj3) == 1

  # rescale: largest dimension → L, origin at (0,0).
  scaled = rescale(proj, 100.0)
  xs = [pt[1] for pt in scaled[1].exterior.points]
  ys = [pt[2] for pt in scaled[1].exterior.points]
  @test minimum(xs) ≈ 0.0   atol = 1e-8
  @test minimum(ys) ≈ 0.0   atol = 1e-8
  @test max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)) ≈ 100.0   atol = 1e-6
end

end # module
