"""
Iberian Peninsula mesh example.

Downloads the Eurostat NUTS-0 (country-level) boundary file, selects ring 1
(mainland) for Spain ("ES") and Portugal ("PT"), unions them into a single
polygon, and produces one combined mesh.

Data source: Eurostat GISCO
  NUTS 2024 – 1:1 Million, EPSG:4326, Level 0
  https://gisco-services.ec.europa.eu/distribution/v2/nuts/
"""

using GeoGmsh
import GeoInterface as GI
import GeometryOps as GO
using LibGEOS          # activates the GeometryOps ↔ GEOS extension
using Downloads

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

nuts0_url  = "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/" *
             "NUTS_RG_01M_2024_4326_LEVL_0.geojson"
nuts0_path = joinpath(data_dir, "NUTS_RG_01M_2024_4326_LEVL_0.geojson")

if !isfile(nuts0_path)
  println("Downloading NUTS-0 boundaries…")
  Downloads.download(nuts0_url, nuts0_path)
  println("  Saved: ", nuts0_path)
end

# ---------------------------------------------------------------------------
# Inspect
# ---------------------------------------------------------------------------

println("\nNUTS-0 Spain + Portugal components (sorted by area):")
comps = list_components(nuts0_path)
comps = filter(row -> row.NUTS_ID ∈ ("ES", "PT"), comps)
sort!(comps, :area, rev = true)
println(comps)

# ---------------------------------------------------------------------------
# Select mainland Spain (ring 1) and mainland Portugal (ring 1)
# ---------------------------------------------------------------------------

geom_col = first(GI.geometrycolumns(comps))

spain_rows    = filter(row -> row.NUTS_ID == "ES" && row.ring == 1, comps)
portugal_rows = filter(row -> row.NUTS_ID == "PT" && row.ring == 1, comps)

spain_geom    = spain_rows[1, geom_col]
portugal_geom = portugal_rows[1, geom_col]

# ---------------------------------------------------------------------------
# Union into one polygon
# ---------------------------------------------------------------------------

println("\nUnioning mainland Spain + Portugal…")
iberia_geom = GO.union(GO.GEOS(), spain_geom, portugal_geom)
println("  Done.  Result trait: ", GI.geomtrait(iberia_geom))

# ---------------------------------------------------------------------------
# 2D geo + mesh
# ---------------------------------------------------------------------------

output = joinpath(data_dir, "iberia")

geoms_to_geo(
  iberia_geom, output;
  target_crs   = "EPSG:3857",
  simplify_tol = 10_000.0,
  bbox_size    = 100.0,
  verbose      = true,
)

geoms_to_msh(
  iberia_geom, output;
  target_crs   = "EPSG:3857",
  simplify_tol = 10_000.0,
  bbox_size    = 100.0,
  mesh_size    = 2.0,
  verbose      = true,
)
