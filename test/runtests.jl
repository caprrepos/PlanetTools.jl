## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

using ColorSchemes: colorschemes
using Documenter: DocMeta, doctest
using Logging: Logging
using Makie: Makie
using PlanetTools
using PlanetTools.Constructors
using PlanetTools.DerivePhysicsVariablesModelE
using PlanetTools.PlanetAerosolUtils
using Test

@testset "PlanetTools" begin
    @testset "Constructors" begin
        include("constructors_tests.jl")
    end

    @testset "DerivePhysicsVariablesModelE" begin
        include("derivephysicsvariablesmodele_tests.jl")
    end

    @testset "PlanetAerosolUtils" begin
        include("planetaerosolutils_tests.jl")
    end

    @testset "Doctests" begin
        setup = quote
            using PlanetTools.Constructors, PlanetTools.DerivePhysicsVariablesModelE
            using PlanetTools.PlanetAerosolUtils
            using ColorSchemes, Makie
        end
        DocMeta.setdocmeta!(PlanetTools, :DocTestSetup, setup; recursive=true)
        doctest(PlanetTools; manual=false)
    end
end
