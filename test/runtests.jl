## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

using ColorSchemes: colorschemes
using Documenter: DocMeta, doctest
using Logging: Logging
using Makie: Makie
using PlanetTools
using PlanetTools.Constructors
using Test

@testset "PlanetTools" begin
    @testset "Constructors" begin
        include("constructors_tests.jl")
    end

    @testset "Doctests" begin
        DocMeta.setdocmeta!(
            PlanetTools,
            :DocTestSetup,
            :(using PlanetTools.Constructors, ColorSchemes, Makie);
            recursive=true,
        )
        doctest(PlanetTools; manual=false)
    end
end
