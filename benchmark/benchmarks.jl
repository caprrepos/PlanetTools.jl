## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-26 to 2026-10-27

# Benchmark suite of PlanetTools, run with AirspeedVelocity.jl (see benchmark/README.md):
#     benchpkg PlanetTools --path=. --rev=main,dirty --exeflags="--threads=4"
# or directly with BenchmarkTools:
#     using BenchmarkTools; include("benchmark/benchmarks.jl"); tune!(SUITE); run(SUITE)
#
# The gridded data have the size of one 3D variable of a ROCKE-3D Mars simulation
# (72 × 46 × 40 grid points, 20 years of monthly means). The PlanetStatistics kernels are
# benchmarked with complete data and with missing data in the layout of real output: grid
# columns below the surface topography are NaN in their lowest levels for all time steps.
# All inputs are generated from fixed seeds, so that every run sees the same data.
#
# The submodules are loaded individually (not with `using PlanetTools`), so that the script
# also runs on older revisions when comparing versions.
using BenchmarkTools
using Dates: Date, Month, Year
using PlanetTools.PlanetAerosolUtils: interpolate_extinction_coeff_log
using PlanetTools.PlanetDustUtils: calc_fractional_dust_masses_bft
using PlanetTools.PlanetStatistics
using PlanetTools.StatisticUtils
using Random: Xoshiro

const SUITE = BenchmarkGroup()

const NLON, NLAT, NLEV, NYEAR = 72, 46, 40, 20
const NT = 12 * NYEAR
const NSIM = 10  # Ensemble size.

# --- Inputs ---

# Data of one simulation and variable, `(Sim, Var, Lon, Lat, Lev, Time)`, with values
# around 200 K. With `missing_frac > 0`, that fraction of grid columns is NaN in its lowest
# `nlev_missing` levels.
function make_data(seed; missing_frac=0.0, nlev_missing=10)
    rng = Xoshiro(seed)
    data = 200.0f0 .+ 20.0f0 .* rand(rng, Float32, 1, 1, NLON, NLAT, NLEV, NT)
    for j in 1:NLAT, i in 1:NLON
        if rand(rng) < missing_frac
            data[1, 1, i, j, 1:nlev_missing, :] .= NaN32
        end
    end
    return data
end

const DATA = Dict(
    "complete data" => make_data(1), "missing data" => make_data(1; missing_frac=0.2)
)
const AXYP = rand(Xoshiro(2), NLON, NLAT) .+ 0.5
const DATES = collect(Date(1, 1):Month(1):Date(NYEAR, 12))
const YEARS = [collect(Date(1):Year(1):Date(NYEAR))]
const T_RANGES = [NT]
const MONTH_WEIGHTS = Float64[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
const MW = repeat(MONTH_WEIGHTS, NYEAR)
const ANNUAL_INDICES = [3 NYEAR]  # Years 3 to 20, after spin-up.
const MONTHLY_INDICES = [collect(m:12:NT) for _ in 1:1, m in 1:12]

nans(dims...) = fill(NaN32, dims...)

# --- PlanetStatistics ---

SUITE["PlanetStatistics"] = BenchmarkGroup()
for (case, data) in DATA
    g = SUITE["PlanetStatistics"][case] = BenchmarkGroup()

    zm, wz = nans(1, 1, 1, NLAT, NLEV, NT), nans(1, 1, 1, NLAT, NLEV, NT)
    gm, wg = nans(1, 1, 1, 1, NLEV, NT), nans(1, 1, 1, 1, NLEV, NT)
    g["calc_planet_means!"] = @benchmarkable calc_planet_means!(
        $data, $zm, $gm, $wz, $wg, $AXYP, $T_RANGES; threshold=80
    )

    # Inputs of the following benchmarks, computed once.
    calc_planet_means!(data, zm, gm, wz, wg, AXYP, T_RANGES; threshold=80)
    ann = nans(1, 1, NLON, NLAT, NLEV, NYEAR)
    zann, gann = nans(1, 1, 1, NLAT, NLEV, NYEAR), nans(1, 1, 1, 1, NLEV, NYEAR)
    wa, wza, wga = similar(ann), similar(zann), similar(gann)
    g["calc_annual_mean_ts!"] = @benchmarkable calc_annual_mean_ts!(
        $data,
        $zm,
        $gm,
        $ann,
        $zann,
        $gann,
        $wa,
        $wza,
        $wga,
        $DATES,
        $T_RANGES,
        $YEARS;
        zm_w=($wz),
        gm_w=($wg),
        month_weights=($MW),
        threshold=80,
    )
    calc_annual_mean_ts!(
        data,
        zm,
        gm,
        ann,
        zann,
        gann,
        wa,
        wza,
        wga,
        DATES,
        T_RANGES,
        YEARS;
        zm_w=wz,
        gm_w=wg,
        month_weights=MW,
        threshold=80,
    )

    moa, soa, wm = (nans(1, 1, NLON, NLAT, NLEV) for _ in 1:3)
    g["calc_multiyear_stats! annual"] = @benchmarkable calc_multiyear_stats!(
        $moa,
        $soa,
        $wm,
        nothing,
        nothing,
        $ann,
        $ANNUAL_INDICES,
        Val(:annual);
        data_w=($wa),
        threshold=80,
    )

    clim, sclim, wclim = (nans(1, 1, NLON, NLAT, NLEV, 12) for _ in 1:3)
    g["calc_multiyear_stats! monthly"] = @benchmarkable calc_multiyear_stats!(
        $clim,
        $sclim,
        $wclim,
        nothing,
        nothing,
        $data,
        $MONTHLY_INDICES,
        Val(:monthly);
        threshold=80,
    )

    calc_multiyear_stats!(
        clim,
        sclim,
        wclim,
        nothing,
        nothing,
        data,
        MONTHLY_INDICES,
        Val(:monthly);
        threshold=80,
    )
    anom = similar(data)
    g["calc_anomalies!"] = @benchmarkable calc_anomalies!(
        $anom, $data, $clim, $MONTHLY_INDICES
    )

    calc_anomalies!(anom, data, clim, MONTHLY_INDICES)
    slopes, pvals = nans(1, 1, NLON, NLAT, NLEV), nans(1, 1, NLON, NLAT, NLEV)
    g["calc_linear_trends_sig! autocorrelation"] = @benchmarkable calc_linear_trends_sig!(
        $slopes, $pvals, $anom, nothing, $DATES; autocorrelation=true
    )
end

# Ensemble of multi-year means, `(Sim, Var, Lon, Lat, Lev)`.
let moa = 200.0f0 .+ 20.0f0 .* rand(Xoshiro(3), Float32, NSIM, 1, NLON, NLAT, NLEV),
    ens_mean = nans(1, NLON, NLAT, NLEV),
    ens_std = nans(1, NLON, NLAT, NLEV)

    SUITE["PlanetStatistics"]["ensemble"] = BenchmarkGroup()
    SUITE["PlanetStatistics"]["ensemble"]["calc_ensemble_stats!"] = @benchmarkable(
        calc_ensemble_stats!($ens_mean, $ens_std, $moa)
    )
end

# --- StatisticUtils ---

# Fields `(Lon, Lat, Lev, Time)` of one simulation; weights `(1, Lat)`.
let rng = Xoshiro(4),
    x = DATA["missing data"][1, 1, :, :, :, :],
    y = x .+ rand(rng, Float32, size(x)...),
    lat = range(-88.0, 88.0; length=NLAT),
    w = permutedims(calc_latitude_weights(collect(lat)))

    g = SUITE["StatisticUtils"] = BenchmarkGroup()
    g["calc_weighted_mean global"] = @benchmarkable calc_weighted_mean($x, $w)
    g["calc_weighted_mean dims=(1, 2)"] = @benchmarkable calc_weighted_mean(
        $x, $w; dims=(1, 2)
    )
    g["calc_weighted_std dims=(1, 2)"] = @benchmarkable calc_weighted_std(
        $x, $w; dims=(1, 2)
    )
    g["calc_weighted_pattern_corr dims=(1, 2)"] = @benchmarkable calc_weighted_pattern_corr(
        $x, $y, $w; dims=(1, 2)
    )
end

# --- PlanetAerosolUtils and PlanetDustUtils ---

# Optical table of 30 spectral bands × 50 radii, interpolated to 200 × 100 points.
let rng = Xoshiro(5),
    edges = exp10.(range(-1, 2; length=31)),  # 0.1 to 100 µm.
    bounds = [(edges[i], edges[i + 1]) for i in 1:30],
    radii = exp10.(range(-2, 1.5; length=50)),
    table = rand(rng, 30, 50) .+ 0.1,
    lambda = exp10.(range(-0.9, 1.9; length=200)),
    r = exp10.(range(-1.9, 1.4; length=100))

    g = SUITE["PlanetAerosolUtils"] = BenchmarkGroup()
    g["interpolate_extinction_coeff_log"] = @benchmarkable interpolate_extinction_coeff_log(
        $bounds, $radii, $table, $lambda, $r
    )
end

let bounds = [0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]  # dust mass bin diameters.
    g = SUITE["PlanetDustUtils"] = BenchmarkGroup()
    g["calc_fractional_dust_masses_bft"] = @benchmarkable calc_fractional_dust_masses_bft(
        $bounds
    )
end
