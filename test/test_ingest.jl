module IngestTests

using GeoGmsh
import GeoInterface as GI
import GeoJSON
using Test

const FIXTURE_SHP = joinpath(@__DIR__, "..", "meshes", "toy", "fixture.shp")

# ---------------------------------------------------------------------------
# Helpers to build raw GI geometries
# ---------------------------------------------------------------------------

# A simple CCW square [0,2]×[0,2]
_square_ccw() = GI.Polygon([GI.LinearRing([(0.,0.),(2.,0.),(2.,2.),(0.,2.),(0.,0.)])])
# The same square wound CW
_square_cw()  = GI.Polygon([GI.LinearRing([(0.,0.),(0.,2.),(2.,2.),(2.,0.),(0.,0.)])])

function run()
  _test_polygon()
  _test_multipolygon()
  _test_feature()
  _test_feature_collection()
  _test_from_dataframe()
  _test_orientation()
  _test_float32_conversion()
end

function _test_polygon()
  geoms = ingest(_square_ccw())
  @test length(geoms) == 1
  g = geoms[1]
  @test g.exterior.closed
  @test npoints(g.exterior) == 4   # closing point stripped
  @test isempty(g.holes)
end

function _test_multipolygon()
  mp    = GI.MultiPolygon([_square_ccw(), _square_cw()])
  geoms = ingest(mp)
  @test length(geoms) == 2
  for g in geoms
    @test npoints(g.exterior) == 4
  end
end

function _test_feature()
  feat  = GI.Feature(_square_ccw(); properties = (name = "test",))
  geoms = ingest(feat)
  @test length(geoms) == 1
end

function _test_feature_collection()
  fc = GI.FeatureCollection([
    GI.Feature(_square_ccw()),
    GI.Feature(_square_cw()),
  ])
  geoms = ingest(fc)
  @test length(geoms) == 2
end

function _test_from_dataframe()
  df = read_geodata(FIXTURE_SHP)
  geoms = ingest(df)
  @test length(geoms) == 3
  for g in geoms
    @test g.exterior.closed
    @test npoints(g.exterior) >= 3
  end
end

function _test_orientation()
  # Regardless of input winding, exterior must be CCW (positive area)
  # and holes must be CW (negative area).
  for poly in [_square_ccw(), _square_cw()]
    g = only(ingest(poly))
    ext_pts = g.exterior.points
    area = GeoGmsh._ring_signed_area(ext_pts)
    @test area > 0   # CCW
  end

  # Polygon with a hole: exterior CCW, hole CW
  outer = GI.LinearRing([(0.,0.),(4.,0.),(4.,4.),(0.,4.),(0.,0.)])
  inner = GI.LinearRing([(1.,1.),(1.,3.),(3.,3.),(3.,1.),(1.,1.)])  # CCW hole (wrong way)
  poly  = GI.Polygon([outer, inner])
  g     = only(ingest(poly))
  @test GeoGmsh._ring_signed_area(g.exterior.points) > 0        # exterior CCW
  @test GeoGmsh._ring_signed_area(g.holes[1].points) < 0        # hole CW
end

function _test_float32_conversion()
  # NaturalEarth returns Float32 — ingest must promote to Float64.
  ring  = GI.LinearRing(NTuple{2,Float32}[(0f0,0f0),(2f0,0f0),(2f0,2f0),(0f0,2f0),(0f0,0f0)])
  poly  = GI.Polygon([ring])
  g     = only(ingest(poly))
  @test eltype(g.exterior.points[1]) == Float64
end

end # module
