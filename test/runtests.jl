## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

using ColorSchemes: colorschemes
using Dates: Dates, Date, Month, Year, daysinmonth, month, year
using Documenter: DocMeta, doctest
using Logging: Logging
using Makie: Makie
using NCDatasets: NCDatasets, NCDataset
using PlanetTools
using PlanetTools.Constructors
using PlanetTools.DerivePhysicsVariables
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

    @testset "DerivePhysicsVariables" begin
        include("derivephysicsvariables_tests.jl")
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

    @testset "Deprecated submodule name" begin
        # DerivePhysicsVariablesModelE is the old name of DerivePhysicsVariables. Binding
        # deprecations print to stderr, not to the logger, so the warning is suppressed here.
        @test Base.isdeprecated(PlanetTools, :DerivePhysicsVariablesModelE)
        @test !(:DerivePhysicsVariablesModelE in names(PlanetTools))
        old = redirect_stderr(devnull) do
            getglobal(PlanetTools, :DerivePhysicsVariablesModelE)
        end
        @test old === PlanetTools.DerivePhysicsVariables
        m = Module()
        redirect_stderr(devnull) do
            Core.eval(m, :(using PlanetTools.DerivePhysicsVariablesModelE))
        end
        @test getglobal(m, :accumulate_variables) === accumulate_variables
    end

    @testset "Deprecated location of compose_palette" begin
        # compose_palette moved from Constructors to PlottingTools after v0.2.1.
        C = PlanetTools.Constructors
        @test Base.isdeprecated(C, :compose_palette)
        old = redirect_stderr(devnull) do
            getglobal(C, :compose_palette)
        end
        @test old === PlanetTools.PlottingTools.compose_palette
        m = Module()
        redirect_stderr(devnull) do
            Core.eval(m, :(using PlanetTools.Constructors: compose_palette))
        end
        @test getglobal(m, :compose_palette) === compose_palette
    end

    @testset "Exports" begin
        # PlanetTools exports every name exported by its submodules.
        submodules = (
            PlanetTools.Constructors,
            PlanetTools.DerivePhysicsVariables,
            PlanetTools.PlanetAerosolUtils,
            PlanetTools.PlanetDustUtils,
            PlanetTools.PlanetStatistics,
            PlanetTools.PlottingTools,
            PlanetTools.ReadWritePlanet,
            PlanetTools.StatisticUtils,
        )
        exported = Set(names(PlanetTools))
        for m in submodules, n in names(m)
            (n === nameof(m) || Base.isdeprecated(m, n)) && continue
            @test n in exported
            @test getfield(PlanetTools, n) === getfield(m, n)
        end
    end

    @testset "CITATION.bib version" begin
        # The citation must name the current version; update CITATION.bib with each release.
        bib = read(joinpath(pkgdir(PlanetTools), "CITATION.bib"), String)
        m = match(r"version\s*=\s*\{v?([^}]+)\}", bib)
        @test m !== nothing
        @test VersionNumber(m.captures[1]) == pkgversion(PlanetTools)
    end

    @testset "Doctests" begin
        setup = quote
            using PlanetTools.Constructors, PlanetTools.DerivePhysicsVariables
            using PlanetTools.PlanetAerosolUtils, PlanetTools.PlanetDustUtils
            using PlanetTools.PlanetStatistics, PlanetTools.PlottingTools
            using PlanetTools.ReadWritePlanet, PlanetTools.StatisticUtils
            using ColorSchemes, Dates, Makie
        end
        DocMeta.setdocmeta!(PlanetTools, :DocTestSetup, setup; recursive=true)
        doctest(PlanetTools; manual=false)
    end
end
