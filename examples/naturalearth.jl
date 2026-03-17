"""
NaturalEarth example — no download required.

Uses NaturalEarth.jl to load country boundaries directly into Julia
(no shapefiles to manage) and produces meshes for France and China.
"""

using GeoGmsh
using NaturalEarth

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

countries = naturalearth("admin_0_countries", 10)

# ---------------------------------------------------------------------------
# France
# ---------------------------------------------------------------------------

println("=== France ===")

geoms_to_geo(
  countries, joinpath(data_dir, "france");
  select       = row -> get(row, :NAME, "") == "France",
  target_crs   = "EPSG:3857",
  simplify_tol = 5_000.0,
  bbox_size    = 100.0,
  verbose      = true,
)

geoms_to_msh(
  countries, joinpath(data_dir, "france");
  select       = row -> get(row, :NAME, "") == "France",
  target_crs   = "EPSG:3857",
  simplify_tol = 5_000.0,
  bbox_size    = 100.0,
  mesh_size    = 2.0,
  verbose      = true,
)

# ---------------------------------------------------------------------------
# China
# ---------------------------------------------------------------------------

println("\n=== China ===")

geoms_to_geo(
  countries, joinpath(data_dir, "china");
  select       = row -> get(row, :NAME, "") == "China",
  target_crs   = "EPSG:3857",
  simplify_tol = 20_000.0,
  bbox_size    = 100.0,
  verbose      = true,
)

geoms_to_msh(
  countries, joinpath(data_dir, "china");
  select       = row -> get(row, :NAME, "") == "China",
  target_crs   = "EPSG:3857",
  simplify_tol = 20_000.0,
  bbox_size    = 100.0,
  mesh_size    = 2.0,
  verbose      = true,
)
