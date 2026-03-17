"""
Gmsh-specific geometry filtering.

Coarsening and refinement are now delegated to the GeometryOps ecosystem:
- Simplification: `GeometryOps.simplify(MinEdgeLength(tol=…), geom)`
- Segmentization: `GeometryOps.segmentize(geom; max_distance=…)`

This file retains only `filter_components`, which removes geometrically
degenerate rings that would cause Gmsh meshing failures.
"""

# ============================================================================
# Public API
# ============================================================================

"""
    filter_components(geoms; min_points = 4) -> Vector{Geometry2D}

Remove geometrically degenerate components:

- Drop any `Geometry2D` whose exterior ring has fewer than `min_points`
  vertices.
- Strip any hole with fewer than `min_points` vertices (the exterior is kept).

The default `min_points = 4` removes 3-point (triangular) rings that survive
simplification and can cause Gmsh meshing failures.
"""
function filter_components(
  geoms      :: Vector{Geometry2D};
  min_points :: Int = 4,
) :: Vector{Geometry2D}
  out = Geometry2D[]
  for g in geoms
    npoints(g.exterior) < min_points && continue
    holes = filter(h -> npoints(h) >= min_points, g.holes)
    push!(out, Geometry2D(g.exterior, holes, g.name))
  end
  return out
end
