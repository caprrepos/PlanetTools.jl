# Benchmarks

`benchmarks.jl` defines the benchmark suite `SUITE` (a BenchmarkTools `BenchmarkGroup`) of
PlanetTools; `plot/plot_benchmarks.jl` draws the results.

## What is measured

| Group | Benchmarks |
|:--|:--|
| `PlanetStatistics / complete data` | `calc_planet_means!`, `calc_annual_mean_ts!`, `calc_multiyear_stats!` (annual and monthly), `calc_anomalies!`, `calc_linear_trends_sig!` with autocorrelation correction |
| `PlanetStatistics / missing data` | the same with missing data |
| `PlanetStatistics / ensemble` | `calc_ensemble_stats!` for 10 simulations |
| `StatisticUtils` | `calc_weighted_mean` (global and per level and time step), `calc_weighted_std`, `calc_weighted_pattern_corr` |
| `PlanetAerosolUtils` | `interpolate_extinction_coeff_log`, 30 bands × 50 radii to 200 × 100 points |
| `PlanetDustUtils` | `calc_fractional_dust_masses_bft` for 8 size bins |

The gridded data have the size of one 3D variable of a ROCKE-3D Mars simulation:
72 × 46 × 40 grid points and 20 years of monthly means (127 MB as `Float32`). In the
"missing data" case, 20 % of the grid columns are `NaN` in their 10 lowest levels for all
time steps, as below the surface topography in real output. All inputs are generated from
fixed seeds.

AirspeedVelocity adds `time_to_load`, the time of `using PlanetTools` (shown as "package
load time").

## Comparing revisions locally

Install AirspeedVelocity once; its commands are placed in `~/.julia/bin`, which must be on
the `PATH`:

```sh
julia -e 'using Pkg; Pkg.add("AirspeedVelocity"); Pkg.build("AirspeedVelocity")'
```

From the repository root, benchmark `main` against the current working copy (`dirty`,
including uncommitted changes):

```sh
benchpkg PlanetTools --path=. --rev=main,dirty --script=benchmark/benchmarks.jl \
    --exeflags="--threads=4" --output-dir=benchmark/results
benchpkgtable PlanetTools --rev=main,dirty --input-dir=benchmark/results --ratio
```

`--script` makes all revisions run the benchmarks of the working copy. Branches, tags, and
commit hashes can be given in `--rev`; `--filter=StatisticUtils` restricts the run to
matching benchmarks. A full run of the suite takes a few minutes per revision.

## Plotting the results

Install the plotting environment once, then plot all results in a directory:

```sh
julia --project=benchmark/plot -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark/plot benchmark/plot/plot_benchmarks.jl \
    --output=benchmark/results.png benchmark/results
```

The figure has three panels with one row per benchmark:

1. **Run time**: median time per call (marker) and the range of the middle 50 % of the
   samples (bar), on a logarithmic axis.
2. **Minimum run time relative to the baseline** (the first revision, or `--baseline=<rev>`):
   below 1 is faster, above 1 slower. The gray band marks ±10 %.
3. **Allocated memory per call.**

The same numbers are written as a Markdown table to `<output>.md`. Options:
`--rev=main,dirty` selects and orders revisions, `--filter=<text>` selects benchmarks, and
`--output` accepts `.png`, `.svg`, and `.pdf`.

To compare thread counts, run `benchpkg` once per thread count into separate directories
and give each result file a label:

```sh
for n in 1 4; do
    benchpkg PlanetTools --path=. --rev=dirty --script=benchmark/benchmarks.jl \
        --exeflags="--threads=$n" --output-dir=benchmark/results/threads$n
done
julia --project=benchmark/plot benchmark/plot/plot_benchmarks.jl \
    --output=benchmark/threads.png \
    "1 thread=benchmark/results/threads1/results_PlanetTools@dirty.json" \
    "4 threads=benchmark/results/threads4/results_PlanetTools@dirty.json"
```

## Directly with BenchmarkTools

BenchmarkTools is not a dependency of PlanetTools; use a separate environment with both:

```julia
using Pkg
Pkg.activate("@planettools-bench")  # A shared environment in ~/.julia/environments.
Pkg.develop(path=".")
Pkg.add("BenchmarkTools")

using BenchmarkTools
include("benchmark/benchmarks.jl")
tune!(SUITE)
results = run(SUITE; verbose=true)
results["PlanetStatistics"]["missing data"]["calc_planet_means!"]
```

## In CI

`.github/workflows/benchmark.yml` benchmarks every pull request to `main` against `main`,
with 2 threads. The comparison table appears in the job summary of the workflow run; the
figure, the table, and the raw results are attached to the run as the artifact
`benchmark-results`. The benchmark script is taken from `main`, so that both revisions run
the same benchmarks; benchmarks added in a pull request take effect after it is merged.

## Reading the results

- **Compare minimum times** between revisions: the minimum is the most robust estimate of
  the run time of the code itself. Medians and spreads include interference from other
  processes and garbage collection.
- **Expect noise on shared machines.** Running the same code twice on a shared 2-core
  cloud machine gave ratios between 0.7 and 1.7 for single benchmarks. On GitHub runners,
  treat differences below about 20 % as noise unless they recur; on a quiet workstation,
  differences of a few percent are meaningful.
- **Allocations are deterministic.** A change in the number of allocations or in the
  allocated memory is always real and often points to a type instability.
- **Keep the thread count fixed** and at most the number of physical cores, and run
  revisions on the same machine in the same session (as `benchpkg` does).
