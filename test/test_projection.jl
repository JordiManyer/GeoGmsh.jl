module ProjectionTests

using GeoGmsh
import GeometryOps as GO
import GeoInterface as GI
using Test

# Small square near lon=10°, lat=20° in WGS84 (degrees).
# A 1° side ≈ 111 km.
const _POLY = GI.Polygon([GI.LinearRing([
  (10.0, 20.0), (11.0, 20.0), (11.0, 21.0), (10.0, 21.0), (10.0, 20.0),
])])

function run()
  # GO.reproject: degrees (EPSG:4326) → metres (EPSG:3857)
  proj = GO.reproject(_POLY; source_crs = "EPSG:4326", target_crs = "EPSG:3857")
  @test !isnothing(proj)
  # At lon≈10°, x ≈ R·deg2rad(10) ≈ 1.11e6 m.
  xs = [GI.x(p) for p in GI.getpoint(GI.getexterior(proj))]
  @test all(x -> 1.0e6 < x < 2.0e6, xs)

  # rescale: works on ingested Geometry2D
  geoms  = ingest(proj)
  scaled = rescale(geoms, 100.0)
  pts = scaled[1].exterior.points
  xs2 = [p[1] for p in pts]
  ys2 = [p[2] for p in pts]
  @test minimum(xs2) ≈ 0.0   atol = 1e-8
  @test minimum(ys2) ≈ 0.0   atol = 1e-8
  @test max(maximum(xs2) - minimum(xs2), maximum(ys2) - minimum(ys2)) ≈ 100.0  atol = 1e-6
end

end # module
