## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

# Tests for `PlottingTools`. Imports are in `runtests.jl`.

newaxis() = Makie.Axis(Makie.Figure()[1, 1])
lineplots(ax) = filter(p -> p isa Makie.Lines, ax.scene.plots)
vertices(p) = Makie.to_value(p[1])
finite_x(p) = [v[1] for v in vertices(p) if isfinite(v[1])]
finite_y(p) = [v[2] for v in vertices(p) if isfinite(v[2])]
plotcolor(p) = Makie.to_value(p.color)

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

@testset "finite_contour!" begin
    x, y = collect(0.0:1.0:4.0), collect(0.0:1.0:3.0)
    plane = [xi for xi in x, yi in y]  # Isoline at level L is the vertical line x = L.

    @testset "geometry" begin
        ax = newaxis()
        res = finite_contour!(ax, x, y, plane; levels=[0.5, 2.5])
        @test res.levels == [0.5, 2.5] && res.nsegments == [1, 1]
        pl = lineplots(ax)
        @test length(pl) == 2
        @test all(≈(0.5), finite_x(pl[1])) && all(≈(2.5), finite_x(pl[2]))
        @test extrema(finite_y(pl[1])) == (0.0, 3.0)
    end

    @testset "missing data leave gaps" begin
        z = copy(plane)
        z[3, 3:4] .= NaN  # x = 2, y ≥ 2.
        ax = newaxis()
        res = finite_contour!(ax, x, y, z; levels=[1.5, 2.5, 3.5])
        pl = lineplots(ax)
        @test length(pl) == 3
        for (p, lev) in zip(pl, res.levels)
            @test all(≈(lev), finite_x(p))  # No spurious vertices.
        end
        # Gap for y > 1 next to the missing cells.
        @test maximum(finite_y(pl[1])) <= 1 && maximum(finite_y(pl[2])) <= 1
        @test maximum(finite_y(pl[3])) == 3  # Unaffected level.
        @test any(v -> isnan(v[1]), vertices(pl[1]))  # Pen-up markers present.

        allnan = fill(NaN, size(plane))
        ax = newaxis()
        res = finite_contour!(ax, x, y, allnan; levels=[1.0])
        @test res.nsegments == [0] && isempty(lineplots(ax))
        @test finite_contour!(newaxis(), x, y, allnan).levels == Float64[]
        @test finite_contour!(newaxis(), x, y, ones(5, 4); levels=3).levels == Float64[]
    end

    @testset "minpoints" begin
        z = copy(plane)
        z[3:end, 3:end] .= NaN  # Level 2.5 keeps a 2-vertex piece at y = 0 and 1.
        n2 = finite_contour!(newaxis(), x, y, z; levels=[2.5]).nsegments
        n9 = finite_contour!(newaxis(), x, y, z; levels=[2.5], minpoints=9).nsegments
        @test n2[1] >= 1 && n9 == [0]
        res = finite_contour!(newaxis(), x, y, plane; levels=[1.5], minpoints=Int32(2))
        @test res.nsegments == [1]
    end

    @testset "level specifications" begin
        z = [xi + 0.0 * yi for xi in x, yi in y]
        @test finite_contour!(newaxis(), x, y, z; levels=3).levels == [1.0, 2.0, 3.0]
        @test finite_contour!(newaxis(), x, y, z; levels=1, zrange=(0, 2)).levels == [1.0]
        @test finite_contour!(newaxis(), x, y, z; levels=1.0).levels == [1.0]
        @test finite_contour!(newaxis(), x, y, z; levels=1 // 2).levels == [0.5]
        @test finite_contour!(newaxis(), x, y, z; levels=[2, 1, 2]).levels == [1.0, 2.0]
        @test finite_contour!(newaxis(), x, y, z, [1.5]; levels=[3.5]).levels == [1.5]
        r = finite_contour!(newaxis(), x, y, z; levels=[10.0, 0.5])
        @test r.levels == [0.5, 10.0] && r.nsegments == [1, 0]  # 10 does not occur.
    end

    @testset "colors" begin
        ax = newaxis()
        finite_contour!(ax, x, y, plane; levels=[2.5, 0.5], color=[:red, :blue])
        pl = lineplots(ax)
        @test plotcolor(pl[1]) == Makie.to_color(:blue)  # Level 0.5.
        @test plotcolor(pl[2]) == Makie.to_color(:red)  # Level 2.5.

        ax = newaxis()
        finite_contour!(ax, x, y, plane; levels=[1, 2, 2], color=[:red, :green, :blue])
        @test plotcolor.(lineplots(ax)) == Makie.to_color.([:red, :green])

        ax = newaxis()
        finite_contour!(ax, x, y, plane; levels=[0.5, 2.5], color=:black, linewidth=3)
        @test all(p -> plotcolor(p) == Makie.to_color(:black), lineplots(ax))
        @test all(p -> Makie.to_value(p.linewidth) == 3, lineplots(ax))

        ax = newaxis()
        finite_contour!(ax, x, y, plane; levels=[1.0, 2.0, 3.0], colormap=:plasma)
        cg = Makie.cgrad(:plasma)
        expected = Makie.to_color.([get(cg, 0.0), get(cg, 0.5), get(cg, 1.0)])
        @test plotcolor.(lineplots(ax)) ≈ expected

        ax = newaxis()
        finite_contour!(ax, x, y, plane; levels=[2.0], colorrange=(0, 4), colormap=:plasma)
        @test plotcolor(only(lineplots(ax))) ≈ Makie.to_color(get(cg, 0.5))
    end

    @testset "invalid input" begin
        ax = newaxis()
        @test_throws DimensionMismatch finite_contour!(ax, x, y[1:3], plane)
        @test_throws ArgumentError finite_contour!(ax, x, y, plane; levels=0)
        @test_throws ArgumentError finite_contour!(ax, x, y, plane; levels=[1.0, NaN])
        @test_throws ArgumentError finite_contour!(
            ax, x, y, plane; levels=[1.0], minpoints=0
        )
        @test_throws ArgumentError finite_contour!(ax, x, y, plane; levels=3, zrange=(2, 1))
        @test_throws ArgumentError finite_contour!(
            ax, x, y, plane; levels=[1.0], colorrange=(0, Inf)
        )
        @test_throws ArgumentError finite_contour!(
            ax, x, y, plane; levels=[1.0, 2.0], color=[:red]
        )
    end
end
