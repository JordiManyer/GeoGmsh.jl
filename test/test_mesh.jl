module MeshTests

using GeoGmsh
using Test

# A small synthetic square (10 km side) so meshing is fast.
const _SIDE  = 10_000.0
const _SQUARE = ShapeGeometry(
  Contour([(0.0,0.0),(0.0,_SIDE),(_SIDE,_SIDE),(_SIDE,0.0)], true),
  Contour[],
)

const TMP_NAME = joinpath(tempdir(), "test_shapefile_to_gmsh_mesh")

function run()
  # Single .msh file, linear elements.
  generate_mesh([_SQUARE], TMP_NAME; mesh_size = 2_000.0)
  @test isfile(TMP_NAME * ".msh")
  rm(TMP_NAME * ".msh")

  # Quadratic elements.
  generate_mesh([_SQUARE], TMP_NAME; mesh_size = 2_000.0, order = 2)
  @test isfile(TMP_NAME * ".msh")
  rm(TMP_NAME * ".msh")

  # Recombined (quad) mesh.
  generate_mesh([_SQUARE], TMP_NAME; mesh_size = 2_000.0, recombine = true)
  @test isfile(TMP_NAME * ".msh")
  rm(TMP_NAME * ".msh")

  # Split-components output.
  two_squares = [_SQUARE,
    ShapeGeometry(
      Contour([(20_000.0,0.0),(20_000.0,_SIDE),(20_000.0+_SIDE,_SIDE),(20_000.0+_SIDE,0.0)], true),
      Contour[],
    )]
  generate_mesh(two_squares, TMP_NAME; mesh_size = 2_000.0, split_components = true)
  @test isdir(TMP_NAME)
  files = readdir(TMP_NAME)
  @test length(files) == 2
  @test all(endswith(".msh"), files)
  rm(TMP_NAME; recursive = true)
end

end # module
