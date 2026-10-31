## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

"""
    PlanetTools.PlanetDustUtils

Submodule of PlanetTools with utilities for mineral dust, currently the size distribution
of emitted dust from brittle fragmentation theory [Kok, 2011]. Load it with
`using PlanetTools.PlanetDustUtils`.
"""
module PlanetDustUtils

using QuadGK: quadgk
using SpecialFunctions: erf

export calc_fractional_dust_masses_bft

"""
    calc_fractional_dust_masses_bft(diameter_bounds; kwargs...) -> Vector{Float64}

Return the fraction of emitted dust mass in each size bin between consecutive
`diameter_bounds`, from the volume size distribution of brittle fragmentation theory
[Kok, 2011, Eq. 6],

    dV/dD ∝ [1 + erf(ln(D / D̄s) / (√2 ln σs))] exp(-(D / λ)³).

`D` is the particle **diameter**, as are the bin bounds of the ModelE sectional dust
scheme; do not pass radii, such as the effective radii used for radiative calculations.
The bin fractions are integrals of this distribution over the bins. Mass fractions equal
volume fractions, i.e. the particle density is assumed to be independent of size.

# Arguments
- `diameter_bounds::AbstractVector{<:Real}`: bin bounds as particle diameters (not radii)
    in µm, strictly increasing, `diameter_bounds[1] ≥ 0`. The last bound may be `Inf`.

# Keywords
- `median_diameter::Real=3.4`: median diameter by volume of the parent soil, `D̄s`, in µm.
- `geometric_sd::Real=3.0`: geometric standard deviation of the parent soil, `σs`, `> 1`.
- `crack_length::Real=12.0`: side crack propagation length, `λ`, in µm.
- `normalize::Bool=true`: if `true`, the fractions sum to 1 over the bins; if `false`, they
    are fractions of the total emitted mass over all sizes.
- `method::Symbol=:integral`: `:integral` integrates numerically; `:midpoint` multiplies
    `dV/dD` at the geometric bin mean by the bin width (earlier implementation, kept to
    reproduce old results; inaccurate for wide bins).

The default parameters are those of Kok [2011] for arid soils.

# Returns
- `Vector{Float64}`: one fraction per bin, `length(diameter_bounds) - 1` values.

# Throws
- `ArgumentError`: if there are fewer than two bounds, the bounds are not strictly
    increasing or are negative, a parameter is out of range, or the distribution is zero
    in all bins.

# References
- Kok, J. F. (2011), A scaling theory for the size distribution of emitted dust aerosols
    suggests climate models underestimate the size of the global dust cycle, *Proc. Natl.
    Acad. Sci. U.S.A.*, 108(3), 1016–1021, doi:10.1073/pnas.1014798108.

# Examples
```jldoctest
julia> f = calc_fractional_dust_masses_bft([0.2, 2, 4, 8, 16]);

julia> round.(f; digits=3)
4-element Vector{Float64}:
 0.045
 0.143
 0.388
 0.425

julia> sum(calc_fractional_dust_masses_bft([0, Inf]; normalize=false)) ≈ 1
true
```
"""
function calc_fractional_dust_masses_bft(
    diameter_bounds::AbstractVector{<:Real};
    median_diameter::Real=3.4,
    geometric_sd::Real=3.0,
    crack_length::Real=12.0,
    normalize::Bool=true,
    method::Symbol=:integral,
)
    d = diameter_bounds
    n = length(d)
    n >= 2 || throw(ArgumentError("need at least two diameter bounds, got $n"))
    if !(d[1] >= 0 && all(i -> d[i] < d[i + 1], 1:(n - 1)))
        throw(ArgumentError("diameter bounds must be non-negative and strictly increasing"))
    end
    median_diameter > 0 || throw(ArgumentError("median_diameter must be positive"))
    geometric_sd > 1 || throw(ArgumentError("geometric_sd must be greater than 1"))
    crack_length > 0 || throw(ArgumentError("crack_length must be positive"))

    dvdd(x) = _bft_volume_density(x, median_diameter, geometric_sd, crack_length)
    if method === :integral
        volume = [quadgk(dvdd, d[i], d[i + 1])[1] for i in 1:(n - 1)]
    elseif method === :midpoint
        if !(d[1] > 0 && isfinite(d[n]))
            throw(ArgumentError("method=:midpoint needs positive, finite diameter bounds"))
        end
        volume = [dvdd(sqrt(d[i] * d[i + 1])) * (d[i + 1] - d[i]) for i in 1:(n - 1)]
    else
        msg = "method must be :integral or :midpoint, got :$method"
        throw(ArgumentError(msg))
    end

    total = normalize ? sum(volume) : quadgk(dvdd, 0, Inf)[1]
    total > 0 || throw(ArgumentError("the size distribution is zero in all bins"))
    return Float64.(volume ./ total)
end

# Unnormalized BFT volume size distribution dV/dD (Kok, 2011, Eq. 6 divided by D).
function _bft_volume_density(d, median_diameter, geometric_sd, crack_length)
    x = log(d / median_diameter) / (sqrt(2) * log(geometric_sd))
    return (1 + erf(x)) * exp(-(d / crack_length)^3)
end

end # module PlanetDustUtils
