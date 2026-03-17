module MeshTests

using GeoGmsh
using Test

# A small synthetic square (10 km side) so meshing is fast.
const _SIDE  = 10_000.0
const _SQUARE = Geometry2D(
  Contour([(0.0,0.0),(0.0,_SIDE),(_SIDE,_SIDE),(_SIDE,0.0)], true),
  Contour[],
)

# Synthetic flat DEM covering the square with padding (constant elevation 100 m).
# transform = [x_origin, x_pixel, x_rot, y_origin, y_rot, y_pixel]
# 5 cols × 6 rows, 5 km pixels, top-left at (-5000, 20000).
const _DEM = DEMRaster(
  fill(100.0, 5, 6),
  [-5_000.0, 5_000.0, 0.0, 20_000.0, 0.0, -5_000.0],
  nothing, nothing,
)
const _SQUARE3D = [lift_to_3d(_SQUARE, _DEM)]

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
    Geometry2D(
      Contour([(20_000.0,0.0),(20_000.0,_SIDE),(20_000.0+_SIDE,_SIDE),(20_000.0+_SIDE,0.0)], true),
      Contour[],
    )]
  generate_mesh(two_squares, TMP_NAME; mesh_size = 2_000.0, split_components = true)
  @test isdir(TMP_NAME)
  files = readdir(TMP_NAME)
  @test length(files) == 2
  @test all(endswith(".msh"), files)
  rm(TMP_NAME; recursive = true)

  # Volume mesh: depth kwarg.
  generate_mesh_volume(_SQUARE3D, _DEM, TMP_NAME; mesh_size = 2_000.0, depth = 500.0)
  @test isfile(TMP_NAME * ".msh")
  content = read(TMP_NAME * ".msh", String)
  @test occursin("\$Elements",      content)
  @test occursin("\$Nodes",         content)
  @test occursin("\$PhysicalNames", content)
  @test occursin("\"Volume\"",      content)
  @test occursin("\"Top\"",         content)
  @test occursin("\"Bottom\"",      content)
  @test occursin("\"Sides\"",       content)
  rm(TMP_NAME * ".msh")

  # Volume mesh: z_bottom kwarg (absolute elevation).
  generate_mesh_volume(_SQUARE3D, _DEM, TMP_NAME; mesh_size = 2_000.0, z_bottom = -200.0)
  @test isfile(TMP_NAME * ".msh")
  rm(TMP_NAME * ".msh")

  # Volume mesh: split_components.
  generate_mesh_volume(_SQUARE3D, _DEM, TMP_NAME; mesh_size = 2_000.0, depth = 500.0,
                        split_components = true)
  @test isdir(TMP_NAME)
  @test length(readdir(TMP_NAME)) == length(_SQUARE3D)
  rm(TMP_NAME; recursive = true)

  # ── Adaptivity types ───────────────────────────────────────────────────────

  # Struct construction and nominal_size.
  sa = SlopeAdaptivity(mesh_size_min = 500.0, mesh_size_max = 2_000.0)
  @test nominal_size(sa) == 2_000.0
  @test sa.percentile == 0.95

  bl = BoundaryLayerAdaptivity(mesh_size = 2_000.0, mesh_size_min = 200.0, width = 3_000.0)
  @test nominal_size(bl) == 2_000.0

  ca = ComposedAdaptivity(sa, bl)
  @test nominal_size(ca) == 2_000.0

  # BoundaryLayerAdaptivity in a 2D mesh (uses curve tags).
  generate_mesh([_SQUARE], TMP_NAME;
    mesh_size = BoundaryLayerAdaptivity(mesh_size = 2_000.0, mesh_size_min = 300.0,
                                         width = 1_500.0))
  @test isfile(TMP_NAME * ".msh")
  rm(TMP_NAME * ".msh")

  # SlopeAdaptivity in a 3D terrain mesh (uses DEM gradient).
  generate_mesh(_SQUARE3D, _DEM, TMP_NAME;
    mesh_size = SlopeAdaptivity(mesh_size_min = 500.0, mesh_size_max = 2_000.0))
  @test isfile(TMP_NAME * ".msh")
  rm(TMP_NAME * ".msh")

  # ComposedAdaptivity in a 3D terrain mesh.
  generate_mesh(_SQUARE3D, _DEM, TMP_NAME;
    mesh_size = ComposedAdaptivity(
      SlopeAdaptivity(mesh_size_min = 500.0, mesh_size_max = 2_000.0),
      BoundaryLayerAdaptivity(mesh_size = 2_000.0, mesh_size_min = 300.0, width = 1_500.0),
    ))
  @test isfile(TMP_NAME * ".msh")
  rm(TMP_NAME * ".msh")
end

end # module
