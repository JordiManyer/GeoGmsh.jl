using GeoGmsh
using Test

include("test_io.jl")
include("test_ingest.jl")
include("test_projection.jl")
include("test_minedgelength.jl")
include("test_edge_ops.jl")
include("test_gmsh.jl")
include("test_mesh.jl")
include("test_pipeline.jl")

@testset "GeoGmsh.jl" begin
  @testset "read_geodata / list_components" begin IOTests.run()            end
  @testset "ingest"                          begin IngestTests.run()        end
  @testset "reproject / rescale"             begin ProjectionTests.run()    end
  @testset "MinEdgeLength"                   begin MinEdgeLengthTests.run() end
  @testset "segmentize / filter_components"  begin EdgeOpsTests.run()       end
  @testset "write_geo"                       begin GmshTests.run()          end
  @testset "generate_mesh"                   begin MeshTests.run()          end
  @testset "geoms_to_geo"                    begin PipelineTests.run()      end
end
