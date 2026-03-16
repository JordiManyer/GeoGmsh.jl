using Documenter
using GeoGmsh

makedocs(
  sitename = "GeoGmsh.jl",
  format   = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  modules   = [GeoGmsh],
  pages     = [
    "Home"           => "index.md",
    "Pipeline guide" => "pipeline.md",
    "API reference"  => "api.md",
  ],
  checkdocs = :exports,
  warnonly  = false,
)

deploydocs(
  repo         = "github.com/JordiManyer/GeoGmsh.jl.git",
  devbranch    = "master",
  push_preview = true,
)
