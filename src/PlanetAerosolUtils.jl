## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

"""
    PlanetTools.PlanetAerosolUtils

Submodule of PlanetTools with aerosol utilities for planetary atmospheres: the
size-integrated effective radius of a binned size distribution and log-log interpolation
of optical lookup tables in wavelength and particle radius. Load it with
`using PlanetTools.PlanetAerosolUtils`.
"""
module PlanetAerosolUtils

using Interpolations: Flat, Gridded, Line, Linear, Throw, extrapolate, interpolate

export calc_size_integrated_effective_radius, interpolate_extinction_coeff_log

"""
    calc_size_integrated_effective_radius(mass, density, radius) -> Real

Return the effective radius of an aerosol population given in size bins:
`sum(V) / sum(V ./ radius)` with bin volumes `V = mass ./ density`.

For spherical particles this equals the area-weighted mean radius
`r_eff = ∫ r³ n(r) dr / ∫ r² n(r) dr` [Hansen and Travis, 1974], evaluated with one
representative radius per bin. That radius should be the bin's own effective radius;
other choices (e.g. the geometric mean) only approximate it.

Pass **radii**, not diameters: size bins defined by diameter bounds (as in ModelE) must
be converted first, e.g. `radius = sqrt.(dlow .* dhigh) ./ 2` for the geometric mean.

Units cancel except for `radius`: the result has the units of `radius`, as long as `mass`
and `density` each use one unit for all bins.

# Arguments
- `mass::AbstractVector{<:Real}`: aerosol mass per bin (e.g. column load or mixing ratio).
- `density`: particle density, one value per bin (`AbstractVector{<:Real}`) or a single
    `Real` for all bins.
- `radius::AbstractVector{<:Real}`: representative (effective) radius of each bin, `> 0`.

# Returns
- `Real`: the effective radius; `NaN` if the total mass is zero.

# Throws
- `ArgumentError`: if the lengths of the vectors differ, they are empty, or a radius or
    density is not positive.

# References
- Hansen, J. E., and L. D. Travis (1974), Light scattering in planetary atmospheres,
    *Space Sci. Rev.*, 16(4), 527–610, doi:10.1007/BF00168069.

# Examples
```jldoctest
julia> r = [0.1, 1.0];  # µm

julia> mass = [1.0, 1.0];  # Equal mass in both bins.

julia> calc_size_integrated_effective_radius(mass, 2650.0, r) ≈ 2 / (1 / 0.1 + 1 / 1.0)
true
```
"""
function calc_size_integrated_effective_radius(
    mass::AbstractVector{<:Real},
    density::Union{Real,AbstractVector{<:Real}},
    radius::AbstractVector{<:Real},
)
    n = length(mass)
    n > 0 || throw(ArgumentError("mass must not be empty"))
    if length(radius) != n || (density isa AbstractVector && length(density) != n)
        throw(ArgumentError("mass, density, and radius must have the same length"))
    end
    all(>(0), radius) || throw(ArgumentError("all radii must be positive"))
    all(>(0), density) || throw(ArgumentError("all densities must be positive"))

    volume = mass ./ density
    return sum(volume) / sum(volume ./ radius)
end

"""
    interpolate_extinction_coeff_log(
        lambda_bounds, orig_radii, ext_matrix, target_lambda, target_radii; kwargs...
    ) -> Matrix

Interpolate a table of extinction coefficients (or another positive optical quantity) from
spectral bands and particle radii to new wavelengths and radii, linearly in
`log(value)` versus `log(wavelength)` and `log(radius)`.

Each band is placed at the geometric mean of its bounds, `sqrt(lower * upper)`. Bands and
radii may be given in any order. Table values `≤ eps_val` are set to `eps_val` before the
logarithm is taken, and results `≤ 1.01 * eps_val` are returned as zero. Near such zeros
the log-space interpolation falls off much faster than a linear one would.

# Arguments
- `lambda_bounds::AbstractVector`: `(lower, upper)` wavelength bounds of the `k` bands,
    `0 < lower ≤ upper`.
- `orig_radii::AbstractVector{<:Real}`: the `l` particle radii of the table, `> 0`.
- `ext_matrix::AbstractMatrix{<:Real}`: `k × l` table; `ext_matrix[i, j]` belongs to band
    `i` and radius `j`.
- `target_lambda::AbstractVector{<:Real}`: wavelengths to interpolate to, `> 0`, in the
    units of `lambda_bounds`.
- `target_radii::AbstractVector{<:Real}`: radii to interpolate to, `> 0`, in the units of
    `orig_radii`.

# Keywords
- `eps_val::Real=1e-12`: floor for table values before the logarithm, in table units.
- `extrapolation::Symbol=:linear`: outside the table, `:linear` continues the log-log
    lines (power laws), `:flat` keeps the edge values, `:error` throws a `BoundsError`.

# Returns
- `Matrix`: `length(target_lambda) × length(target_radii)` interpolated values.

# Throws
- `ArgumentError`: if fewer than two bands or radii are given, sizes do not match, bands or
    radii are repeated, or a wavelength, radius, or bound is not positive.

# Examples
```jldoctest
julia> bounds = [(0.2, 0.8), (0.8, 3.2)];  # Band centers 0.4 and 1.6 µm.

julia> table = [4.0 2.0; 1.0 0.5];  # Rows: bands; columns: radii 1 and 2 µm.

julia> interpolate_extinction_coeff_log(bounds, [1.0, 2.0], table, [0.8], [1.0, 2.0])
1×2 Matrix{Float64}:
 2.0  1.0
```
"""
function interpolate_extinction_coeff_log(
    lambda_bounds::AbstractVector,
    orig_radii::AbstractVector{<:Real},
    ext_matrix::AbstractMatrix{<:Real},
    target_lambda::AbstractVector{<:Real},
    target_radii::AbstractVector{<:Real};
    eps_val::Real=1e-12,
    extrapolation::Symbol=:linear,
)
    nl, nr = length(lambda_bounds), length(orig_radii)
    if nl < 2 || nr < 2
        throw(ArgumentError("need at least two bands and two radii, got $nl and $nr"))
    end
    if size(ext_matrix) != (nl, nr)
        msg = "ext_matrix has size $(size(ext_matrix)), expected ($nl, $nr)"
        throw(ArgumentError(msg))
    end
    if !all(b -> 0 < b[1] <= b[2], lambda_bounds)
        throw(ArgumentError("band bounds must satisfy 0 < lower ≤ upper"))
    end
    for (name, x) in (
        ("orig_radii", orig_radii),
        ("target_lambda", target_lambda),
        ("target_radii", target_radii),
    )
        all(>(0), x) || throw(ArgumentError("all values of $name must be positive"))
    end
    eps_val > 0 || throw(ArgumentError("eps_val must be positive"))
    bc = _boundary_condition(extrapolation)

    # Interpolation knots must be strictly increasing.
    loglambda = [log(sqrt(b[1] * b[2])) for b in lambda_bounds]
    logr = log.(orig_radii)
    il, ir = sortperm(loglambda), sortperm(logr)
    if !(allunique(loglambda) && allunique(logr))
        throw(ArgumentError("band centers and radii must be unique"))
    end
    logext = log.(max.(ext_matrix[il, ir], eps_val))

    itp = extrapolate(interpolate((loglambda[il], logr[ir]), logext, Gridded(Linear())), bc)
    result = [exp(itp(log(tl), log(tr))) for tl in target_lambda, tr in target_radii]
    result[result .<= 1.01 * eps_val] .= 0
    return result
end

function _boundary_condition(extrapolation::Symbol)
    extrapolation === :linear && return Line()
    extrapolation === :flat && return Flat()
    extrapolation === :error && return Throw()
    msg = "extrapolation must be :linear, :flat, or :error, got :$extrapolation"
    return throw(ArgumentError(msg))
end

end # module PlanetAerosolUtils
