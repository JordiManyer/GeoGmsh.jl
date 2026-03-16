using Documenter
using ShapefileToGmsh

makedocs(
  sitename = "ShapefileToGmsh.jl",
  format   = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  modules   = [ShapefileToGmsh],
  pages     = [
    "Home"           => "index.md",
    "Pipeline guide" => "pipeline.md",
    "API reference"  => "api.md",
  ],
  checkdocs = :exports,
  warnonly  = false,
)

deploydocs(
  repo         = "github.com/JordiManyer/ShapefileToGmsh.git",
  devbranch    = "main",
  push_preview = true,
)
