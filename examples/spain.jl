"""
Spain mesh example.

Downloads two Eurostat NUTS boundary files:
  - NUTS-0 (country level) to select the Spanish mainland (ring 1 of the "ES"
    polygon — the largest component once islands are split out).
  - NUTS-2 (region level) to select Catalonia ("ES51", ring 1).

Produces two meshes:
  1. spain     — mainland Spain only.
  2. catalonia — Catalonia only.

Data source: Eurostat GISCO
  NUTS 2024 – 1:1 Million, EPSG:4326
  https://gisco-services.ec.europa.eu/distribution/v2/nuts/
"""

using GeoGmsh
using Downloads

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

base_url = "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/"

nuts0_url  = base_url * "NUTS_RG_01M_2024_4326_LEVL_0.geojson"
nuts0_path = joinpath(data_dir, "NUTS_RG_01M_2024_4326_LEVL_0.geojson")

nuts2_url  = base_url * "NUTS_RG_01M_2024_4326_LEVL_2.geojson"
nuts2_path = joinpath(data_dir, "NUTS_RG_01M_2024_4326_LEVL_2.geojson")

if !isfile(nuts0_path)
  println("Downloading NUTS-0 boundaries…")
  Downloads.download(nuts0_url, nuts0_path)
  println("  Saved: ", nuts0_path)
end

if !isfile(nuts2_path)
  println("Downloading NUTS-2 boundaries (~16 MB)…")
  Downloads.download(nuts2_url, nuts2_path)
  println("  Saved: ", nuts2_path)
end

# ---------------------------------------------------------------------------
# Inspect
# ---------------------------------------------------------------------------

println("\nNUTS-0 Spain components (sorted by area):")
comps = list_components(nuts0_path)
comps = filter(row -> row.NUTS_ID == "ES", comps)
sort!(comps, :area, rev = true)
println(comps)

# ---------------------------------------------------------------------------
# 1. Mainland Spain — ring 1 of the NUTS-0 "ES" polygon
# ---------------------------------------------------------------------------

println("\n=== Spain (mainland) ===")

geoms_to_geo(
  nuts0_path, joinpath(data_dir, "spain");
  select       = row -> row.NUTS_ID == "ES" && row.ring == 1,
  target_crs   = "EPSG:3857",
  simplify_tol = 10_000.0,
  bbox_size    = 100.0,
  verbose      = true,
)

geoms_to_msh(
  nuts0_path, joinpath(data_dir, "spain");
  select       = row -> row.NUTS_ID == "ES" && row.ring == 1,
  target_crs   = "EPSG:3857",
  simplify_tol = 10_000.0,
  bbox_size    = 100.0,
  mesh_size    = 2.0,
  verbose      = true,
)

# ---------------------------------------------------------------------------
# 2. Catalonia (ES51) — ring 1 of the NUTS-2 "ES51" polygon
# ---------------------------------------------------------------------------

println("\n=== Catalonia (ES51) ===")

println("\nNUTS-2 Catalonia components (sorted by area):")
comps2 = list_components(nuts2_path)
comps2 = filter(row -> row.NUTS_ID == "ES51", comps2)
sort!(comps2, :area, rev = true)
println(comps2)

geoms_to_geo(
  nuts2_path, joinpath(data_dir, "catalonia");
  select       = row -> row.NUTS_ID == "ES51" && row.ring == 1,
  target_crs   = "EPSG:3857",
  simplify_tol = 10_000.0,
  bbox_size    = 100.0,
  verbose      = true,
)

geoms_to_msh(
  nuts2_path, joinpath(data_dir, "catalonia");
  select       = row -> row.NUTS_ID == "ES51" && row.ring == 1,
  target_crs   = "EPSG:3857",
  simplify_tol = 10_000.0,
  bbox_size    = 100.0,
  mesh_size    = 2.0,
  verbose      = true,
)
