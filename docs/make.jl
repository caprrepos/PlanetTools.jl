## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-25 to 2026-10-25

# Build the documentation and run the doctests:
#     julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
# Makie is loaded so that the package extension with the plotting methods is available.
using ColorSchemes: ColorSchemes
using Documenter: Documenter, DocMeta, deploydocs, makedocs
using Makie: Makie
using PlanetTools

# Shared with test/runtests.jl, which runs the same doctests in the package tests.
const DOCTEST_SETUP = quote
    using PlanetTools.Constructors, PlanetTools.DerivePhysicsVariables
    using PlanetTools.PlanetAerosolUtils, PlanetTools.PlanetDustUtils
    using PlanetTools.PlanetStatistics, PlanetTools.PlottingTools
    using PlanetTools.ReadWritePlanet, PlanetTools.StatisticUtils
    using ColorSchemes, Dates, Makie
end
DocMeta.setdocmeta!(PlanetTools, :DocTestSetup, DOCTEST_SETUP; recursive=true)

const SUBMODULES = [
    "Constructors",
    "DerivePhysicsVariables",
    "PlanetAerosolUtils",
    "PlanetDustUtils",
    "PlanetStatistics",
    "PlottingTools",
    "ReadWritePlanet",
    "StatisticUtils",
]

makedocs(;
    modules=[PlanetTools],
    authors="Jan P. Perlwitz <jan.p.perlwitz@caprllc.com>",
    sitename="PlanetTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", edit_link="main"
    ),
    pages=["Home" => "index.md", "API" => ["api/$m.md" for m in SUBMODULES]],
)

# Deploys to the gh-pages branch only on GitHub Actions, which sets GITHUB_REPOSITORY.
repo = "github.com/" * get(ENV, "GITHUB_REPOSITORY", "OWNER/PlanetTools.jl") * ".git"
deploydocs(; repo=repo, devbranch="main")
