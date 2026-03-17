"""
Australia mesh example.

Downloads the ABS ASGS 2021 boundary file for the whole continent,
reprojects to Web Mercator, simplifies, and produces both a .geo
script and a split .msh file (one mesh per island component).

Data source: Australian Bureau of Statistics (ABS)
  Australian Statistical Geography Standard (ASGS) Edition 3
  https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3
  File: AUS_2021_AUST_SHP_GDA2020.zip
"""

using GeoGmsh
using Downloads

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

zip_url  = "https://www.abs.gov.au/statistics/standards/" *
           "australian-statistical-geography-standard-asgs-edition-3/" *
           "jul2021-jun2026/access-and-downloads/digital-boundary-files/" *
           "AUS_2021_AUST_SHP_GDA2020.zip"
zip_path = joinpath(data_dir, "AUS_2021_AUST_SHP_GDA2020.zip")

if !isfile(zip_path)
  println("Downloading ABS ASGS boundary file…")
  Downloads.download(zip_url, zip_path)
  println("  Saved: ", zip_path)
end

# ArchGDAL can read directly from a ZIP via the /vsizip/ virtual filesystem.
input  = "/vsizip/$zip_path/AUS_2021_AUST_GDA2020.shp"
output = joinpath(data_dir, "australia")

# ---------------------------------------------------------------------------
# Inspect
# ---------------------------------------------------------------------------

comps = list_components(input)
sort!(comps, :area, rev = true)
println(first(comps, 5))

# ---------------------------------------------------------------------------
# 2D geo + mesh
# ---------------------------------------------------------------------------

# ring == 1 selects only the mainland (largest component); islands are ring ≥ 2.
geoms_to_geo(
  input, output;
  select           = row -> row.AUS_CODE21 == "AUS" && row.ring == 1,
  target_crs       = "EPSG:3857",
  simplify_alg     = MinEdgeLength(tol = 100_000.0),
  bbox_size        = 100.0,
  split_components = true,
  verbose          = true,
)

geoms_to_msh(
  input, output;
  select           = row -> row.AUS_CODE21 == "AUS" && row.ring == 1,
  target_crs       = "EPSG:3857",
  simplify_alg     = MinEdgeLength(tol = 100_000.0),
  bbox_size        = 100.0,
  mesh_size        = 2.0,
  split_components = true,
  verbose          = true,
)
