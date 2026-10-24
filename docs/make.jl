using PlanetTools
using Documenter

DocMeta.setdocmeta!(PlanetTools, :DocTestSetup, :(using PlanetTools); recursive=true)

makedocs(;
    modules=[PlanetTools],
    authors="Jan P Perlwitz <jan.p.perlwitz@caprllc.com>",
    sitename="PlanetTools.jl",
    format=Documenter.HTML(;
        canonical="https://caprrepos.github.io/PlanetTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/caprrepos/PlanetTools.jl", devbranch="main")
