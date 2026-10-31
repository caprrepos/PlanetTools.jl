## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-25 to 2026-10-27

"""
    PlanetTools

Tools for post-processing and analyzing output of climate models for Earth and other
planets. The functions are organized in submodules and all exported by PlanetTools,
so `using PlanetTools` makes them available; to load only some, use a submodule, e.g.
`using PlanetTools.StatisticUtils`. Some functions are specifically designed for
ModelE and ROCKE-3D developed at NASA GISS, others are more generally usable.

# Submodules
- [`Constructors`](@ref PlanetTools.Constructors): ModelE output file names and discrete
    color palettes.
- [`DerivePhysicsVariables`](@ref PlanetTools.DerivePhysicsVariables): sums, differences,
     and quotients of variables, and cumulative dust size classes.
- [`PlanetAerosolUtils`](@ref PlanetTools.PlanetAerosolUtils): effective radius of binned
    size distributions and log-log interpolation of optical lookup tables.
- [`PlanetDustUtils`](@ref PlanetTools.PlanetDustUtils): size distribution of emitted dust
    from brittle fragmentation theory.
- [`PlanetStatistics`](@ref PlanetTools.PlanetStatistics): spatial means, time series,
    multi-year statistics, anomalies, trends, and ensemble statistics.
- [`PlottingTools`](@ref PlanetTools.PlottingTools): isolines of fields with missing data.
- [`ReadWritePlanet`](@ref PlanetTools.ReadWritePlanet): NetCDF input and output.
- [`StatisticUtils`](@ref PlanetTools.StatisticUtils): weighted statistics of gridded
    fields, correlation significance, and latitude weights.

`compose_palette` and `finite_contour!` are defined in a package extension and require
Makie to be loaded (`using Makie`, `using CairoMakie`, …).

# Examples
```jldoctest
julia> using PlanetTools.Constructors

julia> fnames_modele(2010, 12, 2011, 1, "aij", ["E1"])
1×2 Matrix{String}:
 "DEC2010.aijE1.nc"  "JAN2011.aijE1.nc"
```
"""
module PlanetTools

# Include submodules in dependency order
include("Constructors.jl")
include("DerivePhysicsVariables.jl")
include("PlanetAerosolUtils.jl")
include("PlanetDustUtils.jl")
include("PlanetStatistics.jl")
include("PlottingTools.jl")
include("ReadWritePlanet.jl")
include("StatisticUtils.jl")

# Bring the functions of the submodules into the PlanetTools namespace and export them, so
# that `using PlanetTools` provides all of them. Keep the export lists in sync with those of
# the submodules (checked in test/runtests.jl).
using .Constructors
using .DerivePhysicsVariables
using .PlanetAerosolUtils
using .PlanetDustUtils
using .PlanetStatistics
using .PlottingTools
using .ReadWritePlanet
using .StatisticUtils

# Constructors
export fnames_modele, compose_palette
# DerivePhysicsVariables
export accumulate_variables, subtract_variables, divide_variables
export accumulate_dust_masses_modele, accumulate_dust_optics_modele
# PlanetAerosolUtils
export calc_size_integrated_effective_radius, interpolate_extinction_coeff_log
# PlanetDustUtils
export calc_fractional_dust_masses_bft
# PlanetStatistics
export calc_planet_means!, calc_annual_mean_ts!, calc_seasonal_mean_ts!
export calc_multiyear_stats!, calc_anomalies!
export calc_linear_trends!, calc_linear_trends_sig!, generate_significance_mask!
export calc_ensemble_stats!, calc_ensemble_agreement!
# PlottingTools
export finite_contour!
# ReadWritePlanet
export readnetcdf, writenetcdf
# StatisticUtils
export calc_weighted_mean, calc_weighted_std, calc_weighted_pattern_corr
export fisher_z_transform, fisher_z_inv, calc_correlation_significance
export calc_latitude_weights

end  # module PlanetTools
