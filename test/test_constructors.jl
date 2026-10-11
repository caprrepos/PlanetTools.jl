include("../src/Constructors.jl")
using .Constructors, ColorSchemes, Makie, Test

@testset "fnames_modelE" begin
    f = fnames_modelE(1950, 11, 1951, 2, "aij", ["E001", "E002"])
    @test size(f) == (2, 4)
    @test f[1, 1] == "NOV1950.aijE001.nc"
    @test f[2, 4] == "FEB1951.aijE002.nc"
    show(stdout, MIME"text/plain"(), f); println()
    @test fnames_modelE(Int32(850), 1, 850, 1, "taij", split("A B"); ext = "")[2] == "JAN0850.taijB"
    @test size(fnames_modelE(2000, 1, 2009, 12, "aij", ["x"])) == (1, 120)
    @test_throws ArgumentError fnames_modelE(1951, 1, 1950, 1, "aij", ["E001"])
    @test_throws ArgumentError fnames_modelE(1950, 13, 1951, 1, "aij", ["E001"])
    @test_throws ArgumentError fnames_modelE(1950, 1, 1951, 1, "aij", String[])
end

@testset "compose_palette" begin
    br = [-2, -1, -0.5, 0.5, 1, 2]
    p = compose_palette(br; palette1 = reverse(colorschemes[:Blues]), midcolor = :white,
                        palette2 = :Reds, demarcation = 0)
    @test length(p) == 5 && eltype(p) == Makie.RGBAf
    @test p[3] == Makie.to_color(:white)
    @test length(compose_palette(br; palette1 = :Blues, palette2 = :Reds, demarcation = 0)) == 5
    # edge cases that were wrong / crashing before
    @test length(compose_palette([2, 3, 4, 5]; palette1 = :Blues, palette2 = :Reds)) == 3
    @test (@test_logs (:warn,) compose_palette([2, 3, 4, 5]; palette1 = :Blues,
            palette2 = :Reds, midcolor = :white)) |> length == 3
    @test length(compose_palette([0, 1]; palette1 = :Blues)) == 1          # odd-length scheme, 1 color
    @test length(compose_palette([0, 1, 2]; palette1 = :Blues, palette2 = :Reds, midcolor = :white)) == 2
    q = compose_palette([0, 1, 2, 3]; palette1 = :Blues, midcolor = :white)
    @test length(q) == 3 && q[end] == Makie.to_color(:white)
    r = compose_palette([0, 1, 2, 3]; palette2 = "Reds", midcolor = :gray90)
    @test length(r) == 3 && r[1] == Makie.to_color(:gray90)
    @test length(compose_palette(0:10:100; palette1 = :viridis, midcolor = :white)) == 10
    @test_throws ArgumentError compose_palette([1]; palette1 = :Blues)
    @test_throws ArgumentError compose_palette([1, 2, 3])
    @test_throws ArgumentError compose_palette([3, 2, 1]; palette1 = :Blues)
    # usable in Makie
    fig, ax, cf = contourf(1:10, 1:10, rand(10, 10) .* 4 .- 2; levels = br, colormap = p)
    @test cf isa Makie.AbstractPlot
end
