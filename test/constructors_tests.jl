## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

# Tests for `Constructors`. Imports are in `runtests.jl`. The tests of `compose_palette`,
# which moved to `PlottingTools`, are in `plottingtools_tests.jl`.

@testset "fnames_modele" begin
    @testset "names and layout" begin
        f = fnames_modele(1950, 11, 1951, 2, "aij", ["E001", "E002"])
        @test f isa Matrix{String}
        @test size(f) == (2, 4)
        @test f[1, :] == [
            "NOV1950.aijE001.nc",
            "DEC1950.aijE001.nc",
            "JAN1951.aijE001.nc",
            "FEB1951.aijE001.nc",
        ]
        @test f[2, 3] == "JAN1951.aijE002.nc"
    end

    @testset "month count formula" begin
        for (fy, fm, ly, lm) in ((2000, 1, 2000, 1), (1950, 11, 1952, 2), (0, 12, 3, 1))
            nmonths = 12 * (ly - fy) + (lm - fm) + 1
            @test size(fnames_modele(fy, fm, ly, lm, "aij", ["A", "B", "C"])) ==
                (3, nmonths)
        end
    end

    @testset "all month abbreviations" begin
        f = fnames_modele(2000, 1, 2000, 12, "aij", ["X"]; ext="")
        months = split("JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC")
        @test vec(f) == ["$(m)2000.aijX" for m in months]
    end

    @testset "year padding" begin
        @test only(fnames_modele(0, 1, 0, 1, "aij", ["X"])) == "JAN0000.aijX.nc"
        @test only(fnames_modele(850, 1, 850, 1, "aij", ["X"])) == "JAN0850.aijX.nc"
        @test only(fnames_modele(10000, 1, 10000, 1, "aij", ["X"])) == "JAN10000.aijX.nc"
    end

    @testset "extension and ID types" begin
        @test only(fnames_modele(2000, 1, 2000, 1, "taij", ["X"]; ext="")) ==
            "JAN2000.taijX"
        @test only(fnames_modele(2000, 1, 2000, 1, "aij", ["X"]; ext=".txt")) ==
            "JAN2000.aijX.txt"
        ids = split("E001 E002")  # SubStrings
        @test fnames_modele(2000, 1, 2000, 1, "aij", ids) ==
            fnames_modele(2000, 1, 2000, 1, "aij", ["E001", "E002"])
        @test fnames_modele(Int32(2000), Int8(1), 2000, 1, "aij", ["X"]) ==
            ["JAN2000.aijX.nc";;]
    end

    @testset "invalid input" begin
        @test_throws ArgumentError fnames_modele(2000, 0, 2000, 1, "aij", ["X"])
        @test_throws ArgumentError fnames_modele(2000, 13, 2000, 1, "aij", ["X"])
        @test_throws ArgumentError fnames_modele(2000, 1, 2000, 13, "aij", ["X"])
        @test_throws ArgumentError fnames_modele(-5, 1, -5, 2, "aij", ["X"])
        @test_throws ArgumentError fnames_modele(2001, 1, 2000, 1, "aij", ["X"])
        @test_throws ArgumentError fnames_modele(2000, 2, 2000, 1, "aij", ["X"])
        @test_throws ArgumentError fnames_modele(2000, 1, 2000, 1, "aij", String[])
        @test_throws TypeError fnames_modele(2000, 1, 2000, 1, "aij", ["X"]; ext=:nc)
    end
end
