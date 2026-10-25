## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

# Tests for `StatisticUtils`. Imports are in `runtests.jl`.

# Reference implementations on explicit vectors of valid points (Float64).
function ref_stats(x, y, w)
    keep = .!isnan.(x) .& .!isnan.(y) .& (w .> 0)
    xv, yv, wv = Float64.(x[keep]), Float64.(y[keep]), Float64.(w[keep])
    mx, my = sum(wv .* xv) / sum(wv), sum(wv .* yv) / sum(wv)
    cxy = sum(wv .* (xv .- mx) .* (yv .- my))
    r = cxy / sqrt(sum(wv .* (xv .- mx) .^ 2) * sum(wv .* (yv .- my) .^ 2))
    v1, v2 = sum(wv), sum(wv .^ 2)
    sd = sqrt(v1 * sum(wv .* (xv .- mx) .^ 2) / (v1^2 - v2))
    return (; mean=mx, std=sd, corr=r)
end
ref_stats(x, w) = ref_stats(x, zero(x), w)

@testset "weighted mean, std, and pattern correlation" begin
    rng = MersenneTwister(7)
    nlon, nlat, nt = 6, 5, 3
    x = randn(rng, nlon, nlat, nt)
    y = 0.6 .* x .+ randn(rng, nlon, nlat, nt)
    x[rand(rng, size(x)...) .< 0.1] .= NaN
    y[rand(rng, size(y)...) .< 0.1] .= NaN
    wlat = reshape(rand(rng, nlat) .+ 0.1, 1, nlat)  # Broadcast over lon and time.
    wfull = repeat(wlat, nlon, 1, nt)

    @testset "global reduction" begin
        ref = ref_stats(vec(x), vec(y), vec(wfull))
        @test calc_weighted_pattern_corr(x, y, wlat) ≈ ref.corr
        refx = ref_stats(vec(x), vec(wfull))
        @test calc_weighted_mean(x, wlat) ≈ refx.mean
        @test calc_weighted_std(x, wlat) ≈ refx.std
        @test calc_weighted_mean(x, wfull) ≈ calc_weighted_mean(x, wlat)
    end

    @testset "reduction over dims" begin
        r = calc_weighted_pattern_corr(x, y, wlat; dims=(1, 2))
        m = calc_weighted_mean(x, wlat; dims=(1, 2))
        s = calc_weighted_std(x, wlat; dims=(1, 2))
        @test size(r) == size(m) == size(s) == (1, 1, nt)
        for t in 1:nt
            ref = ref_stats(vec(x[:, :, t]), vec(y[:, :, t]), vec(wfull[:, :, t]))
            @test r[1, 1, t] ≈ ref.corr
            refx = ref_stats(vec(x[:, :, t]), vec(wfull[:, :, t]))
            @test m[1, 1, t] ≈ refx.mean && s[1, 1, t] ≈ refx.std
        end
        m1 = calc_weighted_mean(x, wlat; dims=1)
        @test size(m1) == (1, nlat, nt)
        @test m1[1, 2, 3] ≈ ref_stats(x[:, 2, 3], wfull[:, 2, 3]).mean
        @test calc_weighted_mean(x, wlat; dims=(1, 2, 3))[1] ≈ calc_weighted_mean(x, wlat)
    end

    @testset "joint mask for the correlation" begin
        a = collect(1.0:10.0)
        b = 2 .* a .+ 1
        b[8:10] .= NaN
        @test calc_weighted_pattern_corr(a, b, ones(10)) ≈ 1
        @test calc_weighted_pattern_corr(a, -b, 1.0) ≈ -1
    end

    @testset "properties" begin
        @test calc_weighted_pattern_corr(x, x, wlat) ≈ 1
        @test calc_weighted_pattern_corr(x, 3 .* y .+ 7, wlat) ≈
            calc_weighted_pattern_corr(x, y, wlat)
        @test calc_weighted_pattern_corr(x, y, 5 .* wlat) ≈
            calc_weighted_pattern_corr(x, y, wlat)
        @test calc_weighted_std(x .+ 100, wlat) ≈ calc_weighted_std(x, wlat)
        # Equal weights: the usual sample statistics.
        v = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
        @test calc_weighted_mean(v, 1) == 5
        @test calc_weighted_std(v, ones(8)) ≈ sqrt(32 / 7)
    end

    @testset "precision and element types" begin
        x32 = Float32.(x) .+ 1.0f4  # Large offset: accumulation in Float64 matters.
        @test calc_weighted_mean(x32, Float32.(wlat)) isa Float32
        @test calc_weighted_std(x32, Float32.(wlat)) ≈ calc_weighted_std(x, wlat) rtol =
            1e-3
        @test calc_weighted_pattern_corr(Float32.(x), Float32.(y), wlat) isa Float64
        @test calc_weighted_pattern_corr([1, 2, 3, 4], [1, 3, 2, 4], [1, 1, 1, 1]) ≈ 0.8
        @test calc_weighted_mean(big.([1.0, 2.0]), [1, 3]) isa BigFloat
        xm = Union{Missing,Float64}[1.0, missing, 3.0]
        @test calc_weighted_mean(xm, ones(3)) == 2
    end

    @testset "invalid points and degenerate cases" begin
        @test isnan(calc_weighted_mean([1, 2], [0, 0]))
        @test isnan(calc_weighted_mean([NaN, NaN], [1, 1]))
        @test calc_weighted_mean([1.0, 5.0], [1.0, -1.0]) == 1  # Negative weight ignored.
        @test isnan(calc_weighted_std([3.0, NaN], [1, 1]))  # One valid point.
        @test isnan(calc_weighted_pattern_corr([1.0, 1.0, 1.0], [1.0, 2.0, 3.0], 1))
        m = calc_weighted_mean([1.0 NaN; 2.0 NaN], 1; dims=1)
        @test m[1] == 1.5 && isnan(m[2])
    end

    @testset "invalid input" begin
        @test_throws DimensionMismatch calc_weighted_pattern_corr(x, y[:, :, 1:2], wlat)
        @test_throws DimensionMismatch calc_weighted_mean(x, ones(2, nlat))
        @test_throws DimensionMismatch calc_weighted_mean(x, ones(nlon, nlat, nt, 2))
        @test_throws ArgumentError calc_weighted_mean(x, wlat; dims=4)
        @test_throws ArgumentError calc_weighted_mean(x, wlat; dims=(0, 1))
    end
end

@testset "Fisher z and correlation significance" begin
    r = [-0.9, -0.2, 0.0, 0.5, 0.99]
    @test fisher_z_transform(r) ≈ 0.5 .* log.((1 .+ r) ./ (1 .- r))
    @test fisher_z_inv(fisher_z_transform(r)) ≈ r
    @test fisher_z_transform(1.0) == Inf && fisher_z_transform(-1.0) == -Inf
    @test fisher_z_inv(0.0) == 0

    # z-score 1.959963984540054 is the two-sided 5 % point of the standard normal.
    n = 40
    r05 = tanh(1.959963984540054 / sqrt(n - 3))
    p, z = calc_correlation_significance(r05, n)
    @test z ≈ 1.959963984540054 && p ≈ 0.05
    p, z = calc_correlation_significance([0.0, r05, -r05], n)
    @test p ≈ [1.0, 0.05, 0.05] && z ≈ [0, 1.959963984540054, -1.959963984540054]
    p, _ = calc_correlation_significance(0.99, 100)  # Far tail: no cancellation to 0.
    @test 0 < p < 1e-100
    @test calc_correlation_significance(NaN, 10)[1] |> isnan
    @test_throws ArgumentError calc_correlation_significance(0.5, 3)
end

@testset "calc_latitude_weights" begin
    w = calc_latitude_weights(-90:30:90)
    @test sum(w) ≈ 1 && w ≈ reverse(w)
    @test w[1] ≈ 0 atol = 1e-15
    @test w ./ w[4] ≈ cosd.(-90:30:90)
    @test eltype(calc_latitude_weights(Float32[-60, 0, 60])) == Float32
    @test calc_latitude_weights([-60, 0, 60]) ≈ [0.25, 0.5, 0.25]
    @test calc_latitude_weights([95.0]) == [0.0]  # Clipped; no normalization by zero.
end
