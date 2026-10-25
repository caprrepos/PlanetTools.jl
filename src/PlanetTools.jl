## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-25 to 2026-10-25

"""
    PlanetTools

Tools for post-processing and analyzing output of the climate models GISS ModelE and
ROCKE-3D, for Earth and other planets. The functions live in submodules; load the ones you
need, e.g. `using PlanetTools.StatisticUtils`.

# Submodules
- [`Constructors`](@ref PlanetTools.Constructors): ModelE output file names and discrete
    color palettes.
- [`DerivePhysicsVariablesModelE`](@ref PlanetTools.DerivePhysicsVariablesModelE): sums,
    differences, and quotients of variables, and cumulative dust size classes.
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
include("DerivePhysicsVariablesModelE.jl")
include("PlanetAerosolUtils.jl")
include("PlanetDustUtils.jl")
include("PlanetStatistics.jl")
include("PlottingTools.jl")
include("ReadWritePlanet.jl")
include("StatisticUtils.jl")

# Future modules (uncomment as they are developed)
# include("AerosolOptics.jl")
# include("DustEmission.jl")
# include("MarsDynamics.jl")
# include("RadiativeTransfer.jl")
# include("DataIO.jl")
# include("Diagnostics.jl")
# include("Validation.jl")

# Make submodules publicly available at PlanetTools namespace
using .Constructors

end  # module PlanetTools
