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

# ---------------------------------------------------------------------------

"""
    AngleFilter(; tol)

Simplification algorithm that removes any vertex whose interior angle
(∠ a–b–c, in degrees) is below `tol`.

A zig-zag spike has a very acute interior angle at the tip: the two
neighbouring segments nearly reverse direction, so the vectors b→a and b→c
point almost the same way and the angle between them is close to 0°.
Removing such vertices flattens spikes without touching smooth curves or
genuine corners.

A single pass may leave new acute angles at formerly-adjacent points, so
the algorithm iterates until no further points are removed.

# Example

```julia
import GeometryOps as GO
simplified = GO.simplify(AngleFilter(tol = 20.0), polygon)
```
"""
@kwdef struct AngleFilter <: GO.SimplifyAlg
  number :: Union{Int,    Nothing} = nothing
  ratio  :: Union{Float64,Nothing} = nothing
  tol    :: Union{Float64,Nothing} = nothing   # minimum interior angle (degrees)

  function AngleFilter(number, ratio, tol)
    GO._checkargs(number, ratio, tol)
    new(number, ratio, tol)
  end
end

function GO._simplify(alg::AngleFilter, points::Vector, _preserve_endpoint::Bool)
  n = length(points)
  n <= 2 && return points
  cos_tol = cos(deg2rad(alg.tol))   # remove if cos(angle) > cos_tol  (angle < tol)

  # Iterate until stable — removing a spike can expose new acute angles.
  while true
    kept = fill(true, length(points))
    for i in 2:(length(points) - 1)
      a, b, c = points[i-1], points[i], points[i+1]
      bax, bay = a[1] - b[1], a[2] - b[2]
      bcx, bcy = c[1] - b[1], c[2] - b[2]
      len_ba = sqrt(bax^2 + bay^2)
      len_bc = sqrt(bcx^2 + bcy^2)
      (len_ba == 0 || len_bc == 0) && continue
      cos_angle = (bax * bcx + bay * bcy) / (len_ba * len_bc)
      if cos_angle > cos_tol
        kept[i] = false
      end
    end
    all(kept) && break          # stable — no more points removed
    points = points[kept]
    length(points) < 3 && return points
  end

  length(points) < 3 && return points
  return points
end

# ---------------------------------------------------------------------------

"""
    ComposedAlg(first, second)

Applies two simplification algorithms in sequence: `first` is run first,
its output is passed to `second`.

Use the `∘` operator as syntactic sugar (note: `a ∘ b` applies `b` first):

```julia
alg = MinEdgeLength(tol = 5_000.0) ∘ AngleFilter(tol = 20.0)
# equivalent to: ComposedAlg(AngleFilter(tol = 20.0), MinEdgeLength(tol = 5_000.0))
simplified = GO.simplify(alg, polygon)
```
"""
struct ComposedAlg{A <: GO.SimplifyAlg, B <: GO.SimplifyAlg} <: GO.SimplifyAlg
  first  :: A
  second :: B
end

function GO._simplify(alg::ComposedAlg, points::Vector, preserve_endpoint::Bool)
  points = GO._simplify(alg.first,  points, preserve_endpoint)
           GO._simplify(alg.second, points, preserve_endpoint)
end

Base.:∘(a::GO.SimplifyAlg, b::GO.SimplifyAlg) = ComposedAlg(b, a)
