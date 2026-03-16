# Internal helpers for verbose pipeline output.
# Included after projection.jl (uses _global_bbox) and within the module
# (uses Printf, which is in scope via `using Printf`).

# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

# Integer with thousands separators.
function _fmt(n::Integer)
  s = string(abs(n))
  result = Char[]
  for (i, c) in enumerate(reverse(collect(s)))
    i > 1 && (i - 1) % 3 == 0 && push!(result, ',')
    push!(result, c)
  end
  (n < 0 ? "-" : "") * String(reverse(result))
end

# Short label for a CRS argument (source or target).
_crs_label(s::AbstractString) = begin
  m = match(r"""["']([^"']+)["']""", s)
  isnothing(m) ? (length(s) > 50 ? s[1:50] * "…" : s) : m.captures[1]
end
_crs_label(::Proj.Transformation) = "Proj.Transformation"
_crs_label(x) = _crs_label(string(x))
_crs_label(::Nothing) = "—"

# ---------------------------------------------------------------------------
# Geometry summary
# ---------------------------------------------------------------------------

struct _GeomSummary
  n    :: Int
  pts  :: Int
  xmin :: Float64
  xmax :: Float64
  ymin :: Float64
  ymax :: Float64
end

function _geom_summary(geoms::Vector{ShapeGeometry}) :: _GeomSummary
  if isempty(geoms)
    return _GeomSummary(0, 0, 0.0, 0.0, 0.0, 0.0)
  end
  n   = length(geoms)
  pts = sum(npoints(g.exterior) for g in geoms)
  xmin, xmax, ymin, ymax = _global_bbox(geoms)
  return _GeomSummary(n, pts, xmin, xmax, ymin, ymax)
end

function _print_summary(s::_GeomSummary; units = "")
  u = isempty(units) ? "" : "  ($units)"
  println("  Geometries : $(_fmt(s.n))")
  println("  Points     : $(_fmt(s.pts))")
  @printf("  Bbox       : x ∈ [%.4g, %.4g]  y ∈ [%.4g, %.4g]%s\n",
          s.xmin, s.xmax, s.ymin, s.ymax, u)
end

# Just the bbox line.
function _print_bbox(s::_GeomSummary; units = "")
  u = isempty(units) ? "" : "  ($units)"
  @printf("  Bbox       : x ∈ [%.4g, %.4g]  y ∈ [%.4g, %.4g]%s\n",
          s.xmin, s.xmax, s.ymin, s.ymax, u)
end
