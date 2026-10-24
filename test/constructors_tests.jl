## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

# Tests for `fnames_modele` and `compose_palette`. Imports are in `runtests.jl`.

# Colors of the bands between `levels` drawn by `contourf`, as Makie computes them for its
# polygons. Extension bands (values ±Inf, drawn with `lowclip`/`highclip`) are skipped.
function contourf_band_colors(levels; kwargs...)
    x = range(-3, 3; length=400)
    z = [(levels[end] - levels[1]) * (a + 3) / 6 + levels[1] for a in x, _ in x]
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    plt = Makie.contourf!(ax, x, x, z; levels=levels, kwargs...)
    poly = plt.plots[1]
    cmap = Makie.to_value(poly.colormap)
    lo, hi = Makie.to_value(poly.colorrange)
    vals = sort(filter(isfinite, unique(Makie.to_value(poly.color))))
    return [Makie.RGBAf(cmap[(v - lo) / (hi - lo)]) for v in vals]
end

colors_match(a, b) = length(a) == length(b) && all(isapprox.(a, b; atol=1e-3))

const WHITE = Makie.to_color(:white)
const BLUES = colorschemes[:Blues]
const REDS = colorschemes[:Reds]

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

@testset "compose_palette" begin
    b = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]

    @testset "length and type" begin
        for p1 in (nothing, :Blues), p2 in (nothing, :Reds), mid in (nothing, :white)
            p1 === nothing && p2 === nothing && continue
            for breaks in (b, 0:10:100, [0, 1], [1, 2, 3])
                pal = @test_logs min_level = Logging.Error compose_palette(
                    breaks; palette1=p1, midcolor=mid, palette2=p2, demarcation=0
                )
                @test pal isa Vector{Makie.RGBAf}
                @test length(pal) == length(breaks) - 1
            end
        end
    end

    @testset "diverging palette, demarcation between breaks" begin
        pal = compose_palette(
            b; palette1=reverse(BLUES), midcolor=:white, palette2=:Reds, demarcation=0
        )
        @test pal[1] ≈ Makie.to_color(BLUES[end])  # Darkest blue in the lowest class.
        @test pal[2] ≈ Makie.to_color(BLUES[1])
        @test pal[3] == WHITE
        @test pal[4] ≈ Makie.to_color(REDS[1])
        @test pal[5] ≈ Makie.to_color(REDS[end])  # Darkest red in the highest class.
    end

    @testset "demarcation on a break" begin
        pal = compose_palette(
            [-2, -1, 0, 1, 2];
            palette1=:Blues,
            midcolor=:white,
            palette2=:Reds,
            demarcation=0,
        )
        @test pal[3] == WHITE  # Interval [0, 1], just above the demarcation.
        @test pal[4] ≈ Makie.to_color(get(REDS, 0.5))  # Single sample: midpoint.
    end

    @testset "two palettes without midcolor" begin
        pal = compose_palette(b; palette1=:Blues, palette2=:Reds, demarcation=0)
        @test pal[2] ≈ Makie.to_color(BLUES[end])
        @test pal[3] ≈ Makie.to_color(REDS[1])  # Straddling interval goes to palette2.
    end

    @testset "no straddling interval" begin
        kw = (; palette1=:Blues, midcolor=:white, palette2=:Reds)
        pal = @test_logs (:warn, r"no interval straddles") compose_palette(
            [0, 1, 2, 3]; kw..., demarcation=5
        )
        @test pal == compose_palette([0, 1, 2, 3]; palette1=:Blues)
        pal = @test_logs (:warn, r"no interval straddles") compose_palette(
            [0, 1, 2, 3]; kw..., demarcation=-1
        )
        @test pal == compose_palette([0, 1, 2, 3]; palette2=:Reds)
        @test_logs compose_palette(
            [0, 1, 2, 3]; palette1=:Blues, palette2=:Reds, demarcation=5
        )
    end

    @testset "single palette" begin
        pal = compose_palette(0:10:100; palette1=:viridis, midcolor=:white)
        @test pal[end] == WHITE
        @test pal[1:9] == compose_palette(0:10:90; palette1=:viridis)
        pal = compose_palette(0:4; palette2=:Reds, midcolor=:white)
        @test pal[1] == WHITE
        @test pal[2:4] == compose_palette(0:3; palette2=:Reds)
        @test compose_palette([0, 1]; palette1=:Blues, midcolor=:white) == [WHITE]
        @test compose_palette([0, 1]; palette1=:Blues) == [Makie.to_color(get(BLUES, 0.5))]
        # `demarcation` is ignored with one palette.
        @test compose_palette(b; palette1=:Blues, demarcation=0) ==
            compose_palette(b; palette1=:Blues, demarcation=100)
    end

    @testset "palette specification and defaults" begin
        sym = compose_palette(b; palette1=:Blues, palette2=:Reds)
        @test compose_palette(b; palette1="Blues", palette2="Reds") == sym
        @test compose_palette(b; palette1=BLUES, palette2=REDS) == sym
        @test compose_palette(b; palette1=:Blues, palette2=:Reds, demarcation=1.0) == sym
        @test compose_palette(b; palette1=:Blues, palette2=:Reds, demarcation=nothing) ==
            sym
        @test compose_palette(b; palette1=:Blues, midcolor="gray90") ==
            [compose_palette(b[1:(end - 1)]; palette1=:Blues); Makie.to_color("gray90")]
    end

    @testset "invalid input" begin
        @test_throws ArgumentError compose_palette([1]; palette1=:Blues)
        @test_throws ArgumentError compose_palette(Float64[]; palette1=:Blues)
        @test_throws ArgumentError compose_palette([2, 1]; palette1=:Blues)
        @test_throws ArgumentError compose_palette([0, 1, 1, 2]; palette1=:Blues)
        @test_throws ArgumentError compose_palette([0.0, 1.0, NaN]; palette1=:Blues)
        @test_throws ArgumentError compose_palette([0.0, Inf]; palette1=:Blues)
        @test_throws ArgumentError compose_palette([0, 1, 2]; midcolor=:white)
        @test_throws ArgumentError compose_palette([0, 1, 2]; palette1=:NoSuchScheme)
        @test_throws MethodError compose_palette((0, 1, 2); palette1=:Blues)
    end

    # Guards the usage advice in the docstring against changes in Makie's `contourf`.
    @testset "Makie contourf band colors" begin
        pal = compose_palette(
            b; palette1=reverse(BLUES), midcolor=:white, palette2=:Reds, demarcation=0
        )
        @test colors_match(contourf_band_colors(b; colormap=pal), pal)
        @test colors_match(
            contourf_band_colors(
                b; colormap=pal, extendlow=first(pal), extendhigh=last(pal)
            ),
            pal,
        )
        # `:auto` extensions make Makie resample the colormap (documented pitfall).
        @test !colors_match(
            contourf_band_colors(b; colormap=pal, extendlow=:auto, extendhigh=:auto), pal
        )

        uneven = [0.0, 1.0, 2.0, 3.0, 100.0]
        pal = compose_palette(uneven; palette1=:viridis)
        @test colors_match(contourf_band_colors(uneven; colormap=pal), pal)
    end
end
