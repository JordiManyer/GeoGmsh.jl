using ShapefileToGmsh
using Test

include("test_io.jl")
include("test_projection.jl")
include("test_edge_ops.jl")
include("test_gmsh.jl")
include("test_pipeline.jl")

@testset "ShapefileToGmsh.jl" begin
  @testset "read_shapefile"       begin IOTests.run()         end
  @testset "project_to_meters"    begin ProjectionTests.run() end
  @testset "coarsen/refine edges" begin EdgeOpsTests.run()    end
  @testset "write_geo"            begin GmshTests.run()       end
  @testset "shapefile_to_geo"     begin PipelineTests.run()   end
end
