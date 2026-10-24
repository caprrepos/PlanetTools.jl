## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

# Tests for `PlanetStatistics`. Imports are in `runtests.jl`.

# Reference implementations with explicit loops, independent of the module's helpers.
function ref_wmean(x, w)
    keep = .!isnan.(x) .& (w .> 0)
    return sum(w[keep] .* x[keep]) / sum(w[keep])
end

function ref_wmoments(x, w)
    keep = .!isnan.(x) .& (w .> 0)
    xv, wv = x[keep], w[keep]
    mu = sum(wv .* xv) / sum(wv)
    m(k) = sum(wv .* (xv .- mu) .^ k) / sum(wv)
    return mu, sqrt(m(2)), m(3) / m(2)^1.5, m(4) / m(2)^2 - 3
end

nanfill(dims...) = fill(NaN, dims...)

# Random data with about 10 % NaN.
function randdata(rng, dims...)
    x = randn(rng, dims...)
    x[rand(rng, dims...) .< 0.1] .= NaN
    return x
end

const RNG = MersenneTwister(20261024)

@testset "calc_planet_means!" begin
    nsim, nvar, nlon, nlat, nlev, nt = 2, 2, 4, 3, 2, 5
    data = randdata(RNG, nsim, nvar, nlon, nlat, nlev, nt)
    data[1, 1, :, 2, 1, 1] .= NaN  # A fully missing latitude band.
    axyp = repeat([1.0 2.0 1.0], nlon)  # (Lon, Lat)
    t_ranges = [5, 3]
    zm, wz = nanfill(nsim, nvar, 1, nlat, nlev, nt), nanfill(nsim, nvar, 1, nlat, nlev, nt)
    gm, wg = nanfill(nsim, nvar, 1, 1, nlev, nt), nanfill(nsim, nvar, 1, 1, nlev, nt)
    calc_planet_means!(data, zm, gm, wz, wg, axyp, t_ranges; threshold=0)

    for ns in 1:nsim, nv in 1:nvar, k in 1:nlev, t in 1:t_ranges[ns]
        for j in 1:nlat
            x = data[ns, nv, :, j, k, t]
            nvalid = count(!isnan, x)
            @test wz[ns, nv, 1, j, k, t] ≈ nvalid / nlon
            ref = nvalid == 0 ? NaN : ref_wmean(x, ones(nlon))
            @test isequal(zm[ns, nv, 1, j, k, t], ref) || zm[ns, nv, 1, j, k, t] ≈ ref
        end
        x = data[ns, nv, :, :, k, t]
        @test gm[ns, nv, 1, 1, k, t] ≈ ref_wmean(vec(x), vec(axyp))
        @test wg[ns, nv, 1, 1, k, t] ≈ sum(axyp[.!isnan.(x)]) / sum(axyp)  # Area fraction.
    end
    @test all(isnan, zm[2, :, :, :, :, 4:5])  # Beyond t_ranges: untouched.

    # Threshold: complete fields only.
    full = ones(1, 1, nlon, nlat, 1, 1)
    gap = copy(full)
    gap[1, 1, 1, 1, 1, 1] = NaN
    for (x, expected) in ((full, 1.0), (gap, 0.0))
        z, w = nanfill(1, 1, 1, nlat, 1, 1), nanfill(1, 1, 1, nlat, 1, 1)
        g, v = nanfill(1, 1, 1, 1, 1, 1), nanfill(1, 1, 1, 1, 1, 1)
        calc_planet_means!(x, z, g, w, v, axyp, [1])
        @test v[1] == expected
        @test isnan(g[1]) == (expected == 0)
    end
    @test_throws DimensionMismatch calc_planet_means!(
        data, zm, gm, wz, wg, axyp[1:2, :], t_ranges
    )
end

@testset "calc_annual_mean_ts!" begin
    dates = collect(Date(2000, 1):Month(1):Date(2002, 12))
    nt = length(dates)
    days = Float64.(daysinmonth.(dates))
    data = randdata(RNG, 2, 1, 2, 2, nt)
    data_w = rand(RNG, size(data)...)
    data_w[1, 1, 1, 1, 1] = 0.0  # Zero weight counts as missing.
    zm, gm = data[:, :, 1:1, :, :], data[:, :, 1:1, 1:1, :]
    ny = 3
    outs = [nanfill(size(a)[1:(end - 1)]..., ny) for a in (data, zm, gm, data, zm, gm)]
    t_ranges = [nt, 30]  # Second simulation: last year incomplete.
    full_years = [Date(2000):Year(1):Date(2002), Date(2000):Year(1):Date(2002)]
    calc_annual_mean_ts!(
        data,
        zm,
        gm,
        outs...,
        dates,
        t_ranges,
        full_years;
        data_w=data_w,
        month_weights=days,
        threshold=0,
    )
    ann, ann_w = outs[1], outs[4]
    for ns in 1:2, i in 1:2, j in 1:2, y in 1:(ns == 1 ? 3 : 2)
        t = (12 * (y - 1) + 1):(12 * y)
        x, w = data[ns, 1, i, j, t], data_w[ns, 1, i, j, t]
        @test ann[ns, 1, i, j, y] ≈ ref_wmean(x, w .* days[t])  # NaN months excluded.
        @test ann_w[ns, 1, i, j, y] ≈ count(.!isnan.(x) .& (w .> 0)) / 12
    end
    @test all(isnan, ann[2, :, :, :, 3])
    @test outs[3][1, 1, 1, 1, 2] ≈ ref_wmean(gm[1, 1, 1, 1, 13:24], days[13:24])

    # Threshold: one missing month fails 100 %, passes 90 %.
    x = reshape(collect(1.0:12.0), 1, 1, 1, 1, 12)
    x[1, 1, 1, 1, 12] = NaN
    for (thr, mean, avail) in ((100, NaN, 0.0), (90, 6.0, 11 / 12))
        o = [nanfill(1, 1, 1, 1, 1) for _ in 1:6]
        calc_annual_mean_ts!(
            x, x, x, o..., dates[1:12], [12], [Date(2000):Year(1):Date(2000)]; threshold=thr
        )
        @test isequal(o[1][1], mean) && o[4][1] == avail
    end
end

@testset "calc_seasonal_mean_ts!" begin
    dates = collect(Date(2000, 1):Month(1):Date(2001, 12))
    data = randdata(RNG, 1, 1, 2, 2, 24)
    zm, gm = data[:, :, 1:1, :, :], data[:, :, 1:1, 1:1, :]
    # Season 1 (JJA) for 2 years; season 2 (SON) only in the first year.
    ind = zeros(Int, 1, 2, 2)
    ind[1, 1, :] = [6, 18]
    ind[1, 2, 1] = 9
    out = [nanfill(1, 1, 2, size(a)[3:(end - 1)]..., 2) for a in (data, zm, gm)]
    wout = [nanfill(size(o)...) for o in out]
    calc_seasonal_mean_ts!(data, zm, gm, out..., wout..., ind, dates; threshold=0)
    for i in 1:2, j in 1:2
        @test out[1][1, 1, 1, i, j, 2] ≈ ref_wmean(data[1, 1, i, j, 18:20], ones(3))
        @test out[1][1, 1, 2, i, j, 1] ≈ ref_wmean(data[1, 1, i, j, 9:11], ones(3))
        @test wout[1][1, 1, 2, i, j, 1] ≈ count(!isnan, data[1, 1, i, j, 9:11]) / 3
    end
    @test all(isnan, out[1][1, 1, 2, :, :, 2])  # Index 0: no second SON.
end

@testset "calc_multiyear_stats!" begin
    nt = 40
    data = randdata(RNG, 2, 2, 3, nt)
    data_w = rand(RNG, size(data)...)
    weights = rand(RNG, 2, nt)

    @testset "annual" begin
        moa, soa, w, ska, kua = (nanfill(2, 2, 3) for _ in 1:5)
        idx = [1 40; 5 30]
        calc_multiyear_stats!(
            moa,
            soa,
            w,
            ska,
            kua,
            data,
            idx,
            Val(:annual);
            weights=weights,
            data_w=data_w,
            threshold=0,
        )
        for ns in 1:2, v in 1:2, k in 1:3
            t = idx[ns, 1]:idx[ns, 2]
            x = data[ns, v, k, t]
            ref = ref_wmoments(x, data_w[ns, v, k, t] .* weights[ns, t])
            @test [moa[ns, v, k], soa[ns, v, k], ska[ns, v, k], kua[ns, v, k]] ≈
                collect(ref)
            @test w[ns, v, k] ≈ count(!isnan, x) / length(t)
        end
        # Without skewness/kurtosis outputs.
        m2, s2, w2 = nanfill(2, 2, 3), nanfill(2, 2, 3), nanfill(2, 2, 3)
        calc_multiyear_stats!(
            m2,
            s2,
            w2,
            nothing,
            nothing,
            data,
            idx,
            Val(:annual);
            weights=weights,
            data_w=data_w,
            threshold=0,
        )
        @test m2 ≈ moa && s2 ≈ soa
        @test_throws DimensionMismatch calc_multiyear_stats!(
            m2,
            s2,
            w2,
            nothing,
            nothing,
            data,
            idx,
            Val(:annual);
            data_w=data_w[:, :, :, 1:2],
        )
    end

    @testset "monthly" begin
        idx = [collect(m:12:nt) for _ in 1:2, m in 1:12]
        idx[2, 12] = Int[]  # Skipped.
        moa, soa, w = (nanfill(2, 2, 3, 12) for _ in 1:3)
        calc_multiyear_stats!(
            moa, soa, w, nothing, nothing, data, idx, Val(:monthly); threshold=0
        )
        for ns in 1:2, m in 1:11
            x = data[ns, 2, 1, idx[ns, m]]
            mu, sd, _, _ = ref_wmoments(x, ones(length(x)))
            @test moa[ns, 2, 1, m] ≈ mu && soa[ns, 2, 1, m] ≈ sd
        end
        @test all(isnan, moa[2, :, :, 12])
    end

    @testset "seasonal" begin
        sdata = randdata(RNG, 1, 1, 4, 2, 10)  # (Sim, Var, Season, Lat, Year)
        moa, soa, w = (nanfill(1, 1, 4, 2) for _ in 1:3)
        idx = zeros(Int, 1, 4, 2)
        idx[1, 1, :] = [1, 10]
        idx[1, 3, :] = [2, 6]
        calc_multiyear_stats!(
            moa, soa, w, nothing, nothing, sdata, idx, Val(:seasonal); threshold=0
        )
        for (s, t) in ((1, 1:10), (3, 2:6)), j in 1:2
            @test moa[1, 1, s, j] ≈ ref_wmean(sdata[1, 1, s, j, t], ones(length(t)))
        end
        @test all(isnan, moa[1, 1, [2, 4], :])
    end

    @testset "threshold and weights where data are missing" begin
        x = reshape([1.0, 2.0, NaN], 1, 1, 3)
        moa, soa, w = nanfill(1, 1), nanfill(1, 1), nanfill(1, 1)
        calc_multiyear_stats!(
            moa,
            soa,
            w,
            nothing,
            nothing,
            x,
            [1 3],
            Val(:annual);
            data_w=ones(1, 1, 3),
            threshold=0,
        )
        @test moa[1] == 1.5 && w[1] ≈ 2 / 3
        calc_multiyear_stats!(moa, soa, w, nothing, nothing, x, [1 3], Val(:annual))
        @test isnan(moa[1]) && w[1] == 0
    end
end

@testset "calc_anomalies!" begin
    nt = 30
    data = randdata(RNG, 2, 1, 3, nt)
    idx = [collect(m:12:nt) for _ in 1:2, m in 1:12]
    clim = randn(RNG, 2, 1, 3, 12)
    anom = nanfill(size(data)...)
    calc_anomalies!(anom, data, clim, idx)
    for ns in 1:2, k in 1:3, t in 1:nt
        m = mod1(t, 12)
        @test isequal(anom[ns, 1, k, t], data[ns, 1, k, t] - clim[ns, 1, k, m])
    end
    @test_throws DimensionMismatch calc_anomalies!(anom[:, :, :, 1:2], data, clim, idx)
    @test_throws DimensionMismatch calc_anomalies!(anom, data, clim[:, :, 1:2, :], idx)
end

@testset "linear trends" begin
    # Two-sided 5 % and 1 % critical values of Student's t (standard tables).
    pt = PlanetStatistics._two_sided_t_pvalue
    @test pt(2.228138852, 10) ≈ 0.05 rtol = 1e-6
    @test pt(2.085963447, 20) ≈ 0.05 rtol = 1e-6
    @test pt(4.032142984, 5) ≈ 0.01 rtol = 1e-6
    @test pt(0.0, 7) ≈ 1

    @testset "exact recovery and Float64 accuracy" begin
        dates = collect(Date(2000, 1):Month(1):Date(2000, 12))  # One year only.
        x = Float32[0.5 * (m - 1) / 12 + 3 for m in 1:12]
        s = nanfill(1, 1)
        calc_linear_trends!(s, reshape(x, 1, 1, 12), nothing, dates)
        @test s[1] ≈ 0.5 rtol = 1e-5

        dates = collect(Date(1991, 1):Month(1):Date(2020, 12))
        tf = [year(d) + (month(d) - 1) / 12 for d in dates]
        y = 0.02 .* (tf .- 1991) .+ randn(RNG, 360)
        w = rand(RNG, 360)
        X = [ones(360) tf] .* sqrt.(w)
        ref = (X \ (y .* sqrt.(w)))[2]  # Weighted least squares by QR.
        s, p = nanfill(1, 1), nanfill(1, 1)
        calc_linear_trends_sig!(
            s, p, reshape(Float32.(y), 1, 1, 360), reshape(w, 1, 1, 360), dates
        )
        @test s[1] ≈ ref rtol = 1e-5
        calc_linear_trends!(s, reshape(y, 1, 1, 360), reshape(w, 1, 1, 360), dates)
        @test s[1] ≈ ref rtol = 1e-10
    end

    @testset "p-values" begin
        dates = collect(Date(2000, 1):Month(1):Date(2001, 8))
        tf = [year(d) + (month(d) - 1) / 12 for d in dates]
        y = [0.1 * i + (-1)^i for i in 1:20]
        b = (([ones(20) tf]) \ y)[2]
        r = y .- [ones(20) tf] * ([ones(20) tf] \ y)
        se = sqrt(sum(r .^ 2) / 18 / sum((tf .- sum(tf) / 20) .^ 2))
        s, p = nanfill(2, 1), nanfill(2, 1)
        data = permutedims(reshape([y; reverse(y)], 20, 2, 1), (2, 3, 1))
        calc_linear_trends_sig!(s, p, data, nothing, dates)
        @test s[1] ≈ b
        @test p[1] ≈ pt(b / se, 18)
        @test p[1] ≈ 0.01121989 rtol = 1e-5  # Value cross-checked with Distributions.jl.
        @test s[2] ≈ -b && p[2] ≈ p[1]
    end

    @testset "autocorrelation" begin
        dates = collect(Date(1991, 1):Month(1):Date(2010, 12))
        e = zeros(240)
        for i in 2:240
            e[i] = 0.8 * e[i - 1] + randn(RNG)
        end
        x = reshape(0.01 .* (1:240) .+ e, 1, 1, 240)
        s, p, pa = nanfill(1, 1), nanfill(1, 1), nanfill(1, 1)
        calc_linear_trends_sig!(s, p, x, nothing, dates)
        calc_linear_trends_sig!(s, pa, x, nothing, dates; autocorrelation=true)
        @test pa[1] > p[1]
        # Alternating residuals (r₁ < 0): n_e is capped at n, p-value unchanged.
        y = reshape([0.1 * i + (-1)^i for i in 1:20], 1, 1, 20)
        d20 = collect(Date(2000, 1):Month(1):Date(2001, 8))
        calc_linear_trends_sig!(s, p, y, nothing, d20)
        calc_linear_trends_sig!(s, pa, y, nothing, d20; autocorrelation=true)
        @test pa[1] ≈ p[1]
    end

    @testset "missing data and invalid input" begin
        dates = collect(Date(2000, 1):Month(1):Date(2000, 4))
        x = zeros(3, 1, 4)
        x[1, 1, :] = [1.0, NaN, 3.0, 9.0]  # Weight 0 for April: (Jan, 1) and (Mar, 3).
        x[2, 1, :] = [NaN, 4.0, NaN, 2.0]  # (Feb, 4) and (Apr, 2).
        x[3, 1, :] = [NaN, NaN, 5.0, NaN]  # One valid point.
        w = ones(3, 1, 4)
        w[1, 1, 4] = 0
        s, p = nanfill(3, 1), nanfill(3, 1)
        calc_linear_trends!(s, x, w, dates)
        @test s[1] ≈ 12  # 1 per month.
        @test s[2] ≈ -12
        @test isnan(s[3])
        calc_linear_trends_sig!(s, p, x, w, dates)
        @test all(isnan, s) && all(isnan, p)  # Fewer than 3 valid points.
        @test_throws DimensionMismatch calc_linear_trends!(s, x, w, dates[1:3])
        @test_throws DimensionMismatch calc_linear_trends!(s, x, w[:, :, 1:3], dates)
    end
end

@testset "generate_significance_mask!" begin
    mask = fill(-1.0f0, 2, 2)
    generate_significance_mask!(mask, [0.01 0.05; NaN 0.2])
    @test mask == Float32[1 0; 0 0]
    generate_significance_mask!(mask, [0.01 0.05; NaN 0.2]; alpha=0.1)
    @test mask == Float32[1 1; 0 0]
    @test_throws DimensionMismatch generate_significance_mask!(mask, zeros(3))
end

@testset "calc_ensemble_stats!" begin
    moa = randdata(RNG, 6, 2, 3)
    w_moa = rand(RNG, 6, 2, 3)
    w_moa[1, 1, 1] = 0.0
    m, s, sk, ku = (nanfill(2, 3) for _ in 1:4)
    calc_ensemble_stats!(m, s, moa; ens_ska=sk, ens_kua=ku, w_moa=w_moa)
    for v in 1:2, k in 1:3
        ref = ref_wmoments(moa[:, v, k], w_moa[:, v, k])
        @test [m[v, k], s[v, k], sk[v, k], ku[v, k]] ≈ collect(ref)
    end
    # Threshold on Σw / n_sims; weights of missing members do not count.
    x = reshape([1.0, 3.0, NaN, NaN], 4, 1)
    m1, s1 = nanfill(1), nanfill(1)
    calc_ensemble_stats!(m1, s1, x; w_moa=ones(4, 1), threshold=50)
    @test m1[1] == 2.0
    calc_ensemble_stats!(m1, s1, x; w_moa=ones(4, 1), threshold=60)
    @test isnan(m1[1]) && isnan(s1[1])
end

@testset "calc_ensemble_agreement!" begin
    moa = reshape([0.3, 0.1, -0.2, NaN, 0.5, -0.4, 0.0, -0.1], 4, 2)
    w = ones(4, 2)
    w[1, 2] = 0  # Excludes 0.5.
    p, n, a = zeros(Int, 2), zeros(Int, 2), nanfill(2)
    calc_ensemble_agreement!(p, n, a, moa; w_moa=w)
    @test p == [2, 0] && n == [1, 2]
    @test a ≈ [2 / 3, 2 / 3]  # 0.0 counts as valid but neither sign.
    calc_ensemble_agreement!(p, n, a, moa; threshold=80)
    @test p == [0, 1] && n == [0, 2] && isnan(a[1]) && a[2] ≈ 0.5
end
