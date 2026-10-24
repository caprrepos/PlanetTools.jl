## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

# Tests for `PlottingTools`. Imports are in `runtests.jl`.

newaxis() = Makie.Axis(Makie.Figure()[1, 1])
lineplots(ax) = filter(p -> p isa Makie.Lines, ax.scene.plots)
vertices(p) = Makie.to_value(p[1])
finite_x(p) = [v[1] for v in vertices(p) if isfinite(v[1])]
finite_y(p) = [v[2] for v in vertices(p) if isfinite(v[2])]
plotcolor(p) = Makie.to_value(p.color)

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
