"""
    MinEdgeLength(; tol)

Simplification algorithm that removes any vertex whose distance to the
previously *kept* vertex is less than `tol`.

Unlike `RadialDistance` (which measures to the last *visited* point),
`MinEdgeLength` measures to the last *retained* point, giving a strict
guarantee: no edge in the output is shorter than `tol`.

This is useful for Gmsh mesh quality — very short boundary edges can cause
degenerate mesh elements or meshing failures.

# Example

```julia
import GeometryOps as GO
simplified = GO.simplify(MinEdgeLength(tol = 5_000.0), polygon)
```
"""
@kwdef struct MinEdgeLength <: GO.SimplifyAlg
  number :: Union{Int,    Nothing} = nothing
  ratio  :: Union{Float64,Nothing} = nothing
  tol    :: Union{Float64,Nothing} = nothing

  function MinEdgeLength(number, ratio, tol)
    GO._checkargs(number, ratio, tol)
    # Store squared tolerance (same convention as RadialDistance)
    tol = isnothing(tol) ? tol : tol^2
    new(number, ratio, tol)
  end
end

function GO._simplify(alg::MinEdgeLength, points::Vector, _preserve_endpoint::Bool)
  n = length(points)
  n <= 2 && return points
  tol² = alg.tol   # already squared
  kept = fill(false, n)
  kept[1] = true
  kept[n] = true
  prev = 1
  for i in 2:(n - 1)
    if GO._squared_euclid_distance(Float64, points[i], points[prev]) >= tol²
      kept[i] = true
      prev = i
    end
  end
  result = points[kept]
  # GeoInterface requires a LinearRing to have ≥ 3 points.  If the tolerance
  # is so large that the ring collapses, return the original so that ingest +
  # filter_components can discard it cleanly rather than throwing here.
  length(result) < 3 && return points
  return result
end
