# PlanetTools

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caprrepos.github.io/PlanetTools.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caprrepos.github.io/PlanetTools.jl/dev/)
[![Build Status](https://github.com/caprrepos/PlanetTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/caprrepos/PlanetTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/caprrepos/PlanetTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/caprrepos/PlanetTools.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

The package contains helper functions that can be applied in Julia scripts for data
processing and plotting of output from simulations with planetary models. It has been
developed with focus on simulations done with NASA GISS's ModelE and its exoplanet version,
ROCKE-3D, but its usefulness is not necessarily confined to this. The collection of tools
includes reading and writing NetCDF output, deriving new variables, area-weighted and temporal
statistics, trends and their significance, ensemble statistics, aerosol and dust size
distributions, and plotting of fields with missing data.

## Features

The functions are organized in submodules:

| Submodule | Purpose | Main functions |
|:--|:--|:--|
| `Constructors` | ModelE output file names, discrete color palettes | `fnames_modele`, `compose_palette`¹ |
| `DerivePhysicsVariables` | Sums, differences, and quotients of variables; cumulative dust size classes | `accumulate_variables`, `subtract_variables`, `divide_variables`, `accumulate_dust_masses_modele`, `accumulate_dust_optics_modele` |
| `PlanetAerosolUtils` | Effective radius of binned size distributions; log-log interpolation of optical lookup tables | `calc_size_integrated_effective_radius`, `interpolate_extinction_coeff_log` |
| `PlanetDustUtils` | Size distribution of emitted dust from brittle fragmentation theory (Kok, 2011) | `calc_fractional_dust_masses_bft` |
| `PlanetStatistics` | Zonal and global means, annual and seasonal time series, multi-year statistics, anomalies, trends with significance, ensemble statistics | `calc_planet_means!`, `calc_annual_mean_ts!`, `calc_seasonal_mean_ts!`, `calc_multiyear_stats!`, `calc_anomalies!`, `calc_linear_trends!`, `calc_linear_trends_sig!`, `generate_significance_mask!`, `calc_ensemble_stats!`, `calc_ensemble_agreement!` |
| `PlottingTools` | Isolines of fields with missing data, also on GeoMakie axes | `finite_contour!`¹ |
| `ReadWritePlanet` | NetCDF input and output with metadata, missing values, and planetary calendars | `readnetcdf`, `writenetcdf` |
| `StatisticUtils` | Weighted mean, standard deviation, and pattern correlation of gridded fields; Fisher z-transform and correlation significance; latitude weights | `calc_weighted_mean`, `calc_weighted_std`, `calc_weighted_pattern_corr`, `fisher_z_transform`, `fisher_z_inv`, `calc_correlation_significance`, `calc_latitude_weights` |

¹ Requires Makie (see [Installation](#installation)).

All functions are documented with examples in the
[documentation](https://caprrepos.github.io/PlanetTools.jl/stable/) and in the REPL help,
e.g. `?calc_weighted_mean`.

## Installation

PlanetTools requires Julia 1.10 or later and is registered in the Julia General registry:

```julia
using Pkg
Pkg.add("PlanetTools")
```

The plotting functions `compose_palette` and `finite_contour!` live in a
[package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions))
that is loaded automatically once Makie is loaded. To use them, add a Makie backend, and
ColorSchemes if you want to pass color schemes directly:

```julia
Pkg.add(["CairoMakie", "ColorSchemes"])
```

Without Makie, PlanetTools installs no plotting libraries.

## Getting started

`using PlanetTools` makes all functions available. To load only some of them, use the
submodules or import single functions:

```julia
using PlanetTools                    # all functions
using PlanetTools.StatisticUtils     # calc_weighted_mean, calc_latitude_weights, …
using PlanetTools: readnetcdf        # a single function
```

The examples below import the functions they use explicitly, so that it is clear where
each function comes from.

### Global and zonal means of a field

```julia
using PlanetTools.Constructors: fnames_modele
using PlanetTools.ReadWritePlanet: readnetcdf
using PlanetTools.StatisticUtils: calc_latitude_weights, calc_weighted_mean

# Monthly output files "<MON><YYYY>.<gname><run>.nc" of run E1 for Jan. and Feb. 2000.
datadir = "/your/datadir"
files = fnames_modele(2000, 1, 2000, 2, "aij", ["E1"])  # 1×2 Matrix{String}

r = readnetcdf(joinpath(datadir, files[1, 1]), ["tsurf", "lat"])
tsurf = r.variable_data["tsurf"]                     # (lon, lat)
w = calc_latitude_weights(r.variable_data["lat"])'  # (1, lat), cos(lat), sum 1

calc_weighted_mean(tsurf, w)          # area-weighted global mean
calc_weighted_mean(tsurf, w; dims=1)  # zonal mean, size (1, lat)
```

`r.variable_info` holds the attributes and dimensions of each variable,
`r.global_attributes` the global attributes. Missing values are returned as `missing` (or as `NaN` with
`missing_as_nan=true`) and are ignored by the statistics.

### Mass fractions of emitted dust in size bins

```julia
using PlanetTools.PlanetDustUtils: calc_fractional_dust_masses_bft

bounds = [0.1, 2.0, 4.0, 8.0, 16.0]      # bounds as diameters in µm
calc_fractional_dust_masses_bft(bounds)  # ≈ [0.044, 0.143, 0.388, 0.425]
```

### A map with a diverging palette and isolines over missing data

```julia
using CairoMakie
using ColorSchemes: colorschemes
using PlanetTools.Constructors: compose_palette
using PlanetTools.PlottingTools: finite_contour!

lon, lat = -177.5:5.0:177.5, -87.5:5.0:87.5
anom = [6cosd(φ) * sind(λ) for λ in lon, φ in lat]  # (lon, lat), e.g. an anomaly in K
anom[10:20, 5:10] .= NaN                              # missing data

# One color per interval: blues below -1, white for [-1, 1], reds above 1.
breaks = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
colors = compose_palette(
    breaks;
    palette1=reverse(colorschemes[:Blues]),
    midcolor=:white,
    palette2=:Reds,
    demarcation=0,
)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Longitude", ylabel="Latitude")
cf = contourf!(
    ax, lon, lat, anom;
    levels=breaks, colormap=colors, extendlow=first(colors), extendhigh=last(colors),
)
finite_contour!(ax, lon, lat, anom, [-3, 3]; color=:black)  # gaps where data are missing
Colorbar(fig[1, 2], cf)
save("anomaly.png", fig)
```

## Conventions

- **Missing data** are `NaN` in all computations; `readnetcdf` and `writenetcdf` convert
    between `NaN`/`missing` and the NetCDF `_FillValue`.
- **Array layout in `PlanetStatistics`**: `(Sim, Var, Lon, Lat, [Lev], Time)`; reduced
    dimensions are kept as singletons, e.g. zonal means `(Sim, Var, 1, Lat, [Lev], Time)`.
    Results are written into preallocated arrays (functions ending in `!`); initialize them
    with `NaN`.
- **Weights**: temporal weights (e.g. days per month) and availability weights (fraction
    of valid data in `[0, 1]`) are multiplied; a `threshold` sets the minimum percentage of
    valid data for a result.
- **Standard deviations** are sample standard deviations (Bessel's correction generalized
    to weights) unless `corrected=false` is passed.
- **Particle sizes**: dust size bins as **diameters** (in µm), which
    `calc_fractional_dust_masses_bft` expects. Effective radii use **radii**: convert bin
    diameters before calling `calc_size_integrated_effective_radius`, e.g.
    `radius = sqrt.(dlow .* dhigh) ./ 2`.
- **Planetary calendars**: time axes with units such as `"months since …"` and a
    user-defined CF `calendar` (e.g. Mars) are read and written as numbers, without Earth
    date conversion.
- **Multithreading**: `PlanetStatistics` uses all available threads; start Julia with
    `julia --threads=auto`.

## Development

The package is work in progress, and it is being used in ongoing scientific research
projects of the author.

## License

PlanetTools is licensed under the [MIT License](LICENSE). The licenses of all dependencies,
including the binary libraries they download, are listed in
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md), generated by
[`scripts/license_check.jl`](scripts/license_check.jl); the assessment is in
[`docs/LICENSE_REVIEW.md`](docs/LICENSE_REVIEW.md). Without Makie, the only copyleft
dependency is the LGPL library libiconv (via NetCDF). The optional Makie extension brings
GPL-licensed and non-redistributable binaries (e.g. FFmpeg); this matters only if you
redistribute binaries, e.g. in a sysimage or container.

## Citing

If you use PlanetTools in your work, please cite it with the entry from
[`CITATION.bib`](CITATION.bib) and for the applied methods the
original scientific references.

## References

Hansen, J. E., and L. D. Travis (1974), Light scattering in planetary atmospheres,
*Space Sci. Rev.*, 16, 527–610, doi:[10.1007/BF00168069](https://doi.org/10.1007/BF00168069).

Kok, J. F. (2011), A scaling theory for the size distribution of emitted dust aerosols
suggests climate models underestimate the size of the global dust cycle, *Proc. Natl. Acad.
Sci. USA*, 108(3), 1016–1021,
doi:[10.1073/pnas.1014798108](https://doi.org/10.1073/pnas.1014798108).

Santer, B. D., T. M. L. Wigley, J. S. Boyle, D. J. Gaffen, J. J. Hnilo, D. Nychka,
D. E. Parker, and K. E. Taylor (2000), Statistical significance of trends and trend
differences in layer-average atmospheric temperature time series, *J. Geophys. Res.*,
105(D6), 7337–7356, doi:[10.1029/1999JD901105](https://doi.org/10.1029/1999JD901105).

## Author

Jan P. Perlwitz, CAPR LLC (Climate, Aerosol, and Pollution Research), science collaborator
at NASA GISS.

Parts of the code, tests, and documentation were developed with AI assistance (Claude,
Anthropic); the files concerned say so in their headers.
