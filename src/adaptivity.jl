"""
Edge resolution operations: coarsening and refinement.

## Coarsening

Remove points that form edges shorter than a minimum length threshold.
Two strategies are available:

- **Single-pass** (`:single`): one left-to-right sweep. Fast; may leave
  residual short edges if removing a point creates a new short edge with
  its new neighbour.
- **Iterative** (`:iterative`): repeats single-pass sweeps until the point
  count stabilises (no more removals). Guarantees no edge shorter than
  `min_edge_length` remains.

## Refinement

Insert equispaced intermediate points along edges that exceed a maximum
length. Given an edge of length L > `max_edge_length`, it is split into
`n = ceil(L / max_edge_length)` sub-edges of equal length.
"""

# ============================================================================
# Public API
# ============================================================================

"""
    coarsen_edges(geoms, min_edge_length; strategy=:iterative)
    -> Vector{ShapeGeometry}

Remove points whose distance to the previous kept point is less than
`min_edge_length` (in the same units as the coordinates — use metres after
reprojection).

`strategy` is `:single` (one sweep) or `:iterative` (sweep until stable).
"""
function coarsen_edges(
  geoms           :: Vector{ShapeGeometry},
  min_edge_length :: Real;
  strategy        :: Symbol = :iterative,
) :: Vector{ShapeGeometry}
  len = Float64(min_edge_length)
  return [_coarsen_geom(g, len, strategy) for g in geoms]
end

"""
    refine_edges(geoms, max_edge_length) -> Vector{ShapeGeometry}

Subdivide edges longer than `max_edge_length` by inserting equispaced
intermediate points along straight-line segments.
"""
function refine_edges(
  geoms           :: Vector{ShapeGeometry},
  max_edge_length :: Real,
) :: Vector{ShapeGeometry}
  len = Float64(max_edge_length)
  return [_refine_geom(g, len) for g in geoms]
end

"""
    filter_components(geoms; min_points = 4) -> Vector{ShapeGeometry}

Remove geometrically degenerate components:

- Drop any `ShapeGeometry` whose exterior ring has fewer than `min_points`
  vertices.
- Strip any hole with fewer than `min_points` vertices (the exterior is kept).

The default `min_points = 4` removes 3-point (triangular) rings that survive
coarsening and can cause Gmsh meshing failures.
"""
function filter_components(
  geoms      :: Vector{ShapeGeometry};
  min_points :: Int = 4,
) :: Vector{ShapeGeometry}
  out = ShapeGeometry[]
  for g in geoms
    npoints(g.exterior) < min_points && continue
    holes = filter(h -> npoints(h) >= min_points, g.holes)
    push!(out, ShapeGeometry(g.exterior, holes, g.name))
  end
  return out
end

# Single-geometry overloads for convenience.
coarsen_edges(g::ShapeGeometry, len::Real; kwargs...) =
  only(coarsen_edges([g], len; kwargs...))

refine_edges(g::ShapeGeometry, len::Real) =
  only(refine_edges([g], len))

# ============================================================================
# Coarsening internals
# ============================================================================

function _coarsen_geom(g::ShapeGeometry, len::Float64, strategy::Symbol)
  ShapeGeometry(
    _coarsen_contour(g.exterior, len, strategy),
    [_coarsen_contour(h, len, strategy) for h in g.holes],
    g.name,
  )
end

function _coarsen_contour(c::Contour, len::Float64, strategy::Symbol)
  if strategy == :single
    return _coarsen_single(c, len)
  elseif strategy == :iterative
    return _coarsen_iterative(c, len)
  else
    throw(ArgumentError(
      "Unknown coarsening strategy: $strategy. Use :single or :iterative."
    ))
  end
end

# Single left-to-right sweep: keep a point only if it is at least `len`
# away from the last kept point.
function _coarsen_single(c::Contour, len::Float64)
  pts = c.points
  n   = length(pts)

  # Need at least 3 points for a closed ring, or 2 for an open polyline.
  min_pts = c.closed ? 3 : 2
  n <= min_pts && return c

  kept = NTuple{2,Float64}[pts[1]]
  for i in 2:n
    if _dist(kept[end], pts[i]) >= len
      push!(kept, pts[i])
    end
  end

  # For closed contours: check that the last kept point is not too close
  # to the first one (wrapping edge).
  if c.closed && length(kept) > min_pts
    while length(kept) > min_pts && _dist(kept[end], kept[1]) < len
      pop!(kept)
    end
  end

  # Guard: never reduce below minimum viable ring.
  if length(kept) < min_pts
    return c
  end

  return Contour(kept, c.closed)
end

# Repeat single-pass until point count stabilises.
function _coarsen_iterative(c::Contour, len::Float64)
  current = c
  while true
    next = _coarsen_single(current, len)
    npoints(next) == npoints(current) && break
    current = next
  end
  return current
end

# ============================================================================
# Refinement internals
# ============================================================================

function _refine_geom(g::ShapeGeometry, len::Float64)
  ShapeGeometry(
    _refine_contour(g.exterior, len),
    [_refine_contour(h, len) for h in g.holes],
    g.name,
  )
end

function _refine_contour(c::Contour, max_len::Float64)
  pts   = c.points
  n     = length(pts)
  edges = c.closed ? n : n - 1

  result = NTuple{2,Float64}[]
  sizehint!(result, n)

  for i in 1:edges
    p1 = pts[i]
    p2 = pts[mod1(i + 1, n)]
    push!(result, p1)
    d = _dist(p1, p2)
    if d > max_len
      nsub = ceil(Int, d / max_len)
      for k in 1:(nsub - 1)
        t = k / nsub
        push!(result, (
          p1[1] + t * (p2[1] - p1[1]),
          p1[2] + t * (p2[2] - p1[2]),
        ))
      end
    end
  end

  # For open contours add the final point (closed contours wrap back to [1]).
  if !c.closed
    push!(result, pts[end])
  end

  return Contour(result, c.closed)
end

# ============================================================================
# Utility
# ============================================================================

@inline _dist(a::NTuple{2,Float64}, b::NTuple{2,Float64}) =
  hypot(b[1] - a[1], b[2] - a[2])
