## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-26 to 2026-10-26

# Benchmarks of the PlanetStatistics kernels at the size of one 3D variable of a ROCKE-3D
# Mars simulation (72 × 46 × 40, 20 years of monthly means).
# Run with: using PkgBenchmark; benchmarkpkg("PlanetTools")
using BenchmarkTools
using Dates: Date, Month, Year
using PlanetTools.PlanetStatistics

SUITE = BenchmarkGroup()

const NLON, NLAT, NLEV, NT = 72, 46, 40, 240
const DATA = rand(Float32, 1, 1, NLON, NLAT, NLEV, NT)
const AXYP = rand(NLON, NLAT)
const DATES = collect(Date(1, 1):Month(1):Date(20, 12))
const MW = rand(NT)

zm, wz = fill(NaN32, 1, 1, 1, NLAT, NLEV, NT), fill(NaN32, 1, 1, 1, NLAT, NLEV, NT)
gm, wg = fill(NaN32, 1, 1, 1, 1, NLEV, NT), fill(NaN32, 1, 1, 1, 1, NLEV, NT)
SUITE["calc_planet_means!"] = @benchmarkable calc_planet_means!(
    $DATA, $zm, $gm, $wz, $wg, $AXYP, [NT]; threshold=80
)

ann = fill(NaN32, 1, 1, NLON, NLAT, NLEV, 20)
zann, gann = fill(NaN32, 1, 1, 1, NLAT, NLEV, 20), fill(NaN32, 1, 1, 1, 1, NLEV, 20)
wa, wza, wga = similar(ann), similar(zann), similar(gann)
calc_planet_means!(DATA, zm, gm, wz, wg, AXYP, [NT]; threshold=80)
SUITE["calc_annual_mean_ts!"] = @benchmarkable calc_annual_mean_ts!(
    $DATA,
    $zm,
    $gm,
    $ann,
    $zann,
    $gann,
    $wa,
    $wza,
    $wga,
    $DATES,
    [NT],
    [Date(1):Year(1):Date(20)];
    zm_w=($wz),
    gm_w=($wg),
    month_weights=($MW),
    threshold=80,
)

moa, soa, wm = (fill(NaN32, 1, 1, NLON, NLAT, NLEV) for _ in 1:3)
SUITE["calc_multiyear_stats! annual"] = @benchmarkable calc_multiyear_stats!(
    $moa,
    $soa,
    $wm,
    nothing,
    nothing,
    $ann,
    [3 20],
    Val(:annual);
    data_w=($wa),
    threshold=80,
)
