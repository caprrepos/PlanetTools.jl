## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

using ColorSchemes: colorschemes
using Dates: Date, Month, Year, daysinmonth, month, year
using Documenter: DocMeta, doctest
using Logging: Logging
using Makie: Makie
using NCDatasets: NCDatasets, NCDataset
using PlanetTools
using PlanetTools.Constructors
using PlanetTools.DerivePhysicsVariablesModelE
using PlanetTools.PlanetAerosolUtils
using PlanetTools.PlanetDustUtils
using PlanetTools.PlanetStatistics
using PlanetTools.PlanetStatistics: PlanetStatistics
using PlanetTools.PlottingTools
using PlanetTools.ReadWritePlanet
using PlanetTools.StatisticUtils
using Random: MersenneTwister
using SpecialFunctions: erf
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

    @testset "PlanetDustUtils" begin
        include("planetdustutils_tests.jl")
    end

    @testset "PlanetStatistics" begin
        include("planetstatistics_tests.jl")
    end

    @testset "PlottingTools" begin
        include("plottingtools_tests.jl")
    end

    @testset "ReadWritePlanet" begin
        include("readwriteplanet_tests.jl")
    end

    @testset "StatisticUtils" begin
        include("statisticutils_tests.jl")
    end

    @testset "Doctests" begin
        setup = quote
            using PlanetTools.Constructors, PlanetTools.DerivePhysicsVariablesModelE
            using PlanetTools.PlanetAerosolUtils, PlanetTools.PlanetDustUtils
            using PlanetTools.PlanetStatistics, PlanetTools.PlottingTools
            using PlanetTools.ReadWritePlanet, PlanetTools.StatisticUtils
            using ColorSchemes, Dates, Makie
        end
        DocMeta.setdocmeta!(PlanetTools, :DocTestSetup, setup; recursive=true)
        doctest(PlanetTools; manual=false)
    end
end
