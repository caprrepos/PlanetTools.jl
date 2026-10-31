## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

"""
    PlanetTools.StatisticUtils

Submodule of PlanetTools with weighted statistics of gridded fields (mean, standard
deviation, pattern correlation), Fisher z-transform and correlation significance, and
latitude weights. Load it with `using PlanetTools.StatisticUtils`.

The weighted functions take arrays of any rank and reduce over `dims` (`:` for all
dimensions, an integer, or a tuple). Weights must be broadcast-compatible with the data,
e.g. `(1, nlat)` latitude weights for `(nlon, nlat, …)` data. Points with `NaN` or `missing`
data or with weights `≤ 0` are ignored. Sums are accumulated in at least `Float64`.
"""
module StatisticUtils

using SpecialFunctions: erfc

export calc_weighted_mean, calc_weighted_std, calc_weighted_pattern_corr
export fisher_z_transform, fisher_z_inv, calc_correlation_significance
export calc_latitude_weights

"""
    calc_weighted_mean(vals, weights; dims=:)

Return the weighted mean `Σ w x / Σ w` of `vals` over `dims`, ignoring invalid points.

# Arguments
- `vals::AbstractArray{<:Union{Real,Missing}}`: data of any rank.
- `weights`: weights broadcast-compatible with `vals` (an array or a `Real`).

# Keywords
- `dims=:`: dimensions to reduce; `:` returns a scalar.

# Returns
- The mean (an array with reduced dimensions of size 1, or a scalar for `dims=:`) in the
    floating-point type of the inputs; `NaN` where no point is valid.

See also [`calc_weighted_std`](@ref), [`calc_weighted_pattern_corr`](@ref).

# Examples
```jldoctest
julia> x = [1.0 2.0; NaN 4.0];  # (lon, lat)

julia> w = [1.0 3.0];  # Latitude weights, shape (1, lat).

julia> calc_weighted_mean(x, w) ≈ (1 * 1 + 2 * 3 + 4 * 3) / 7
true

julia> calc_weighted_mean(x, w; dims=1)
1×2 Matrix{Float64}:
 1.0  3.0
```
"""
function calc_weighted_mean(vals::AbstractArray, weights; dims=:)
    keep, w = _prepare(vals, nothing, weights, dims)
    s = _first_moments(vals, nothing, w, keep)
    RT = _rtype(vals, nothing, w)
    return _finish(dims, RT.(_ratio.(s.wx, s.w)))
end

"""
    calc_weighted_std(vals, weights; dims=:)

Return the weighted standard deviation of `vals` over `dims` for reliability weights,
`sqrt(V₁ Σ w (x - μ)² / (V₁² - V₂))` with `V₁ = Σ w` and `V₂ = Σ w²`, ignoring invalid
points. For equal weights this is the usual sample standard deviation (`n - 1`).

# Arguments
- `vals::AbstractArray{<:Union{Real,Missing}}`: data of any rank.
- `weights`: weights broadcast-compatible with `vals` (an array or a `Real`).

# Keywords
- `dims=:`: dimensions to reduce; `:` returns a scalar.

# Returns
- The standard deviation, shaped as for [`calc_weighted_mean`](@ref); `NaN` where fewer
    than two points are valid (`V₁² ≤ V₂`).

# Examples
```jldoctest
julia> calc_weighted_std([1.0, 2.0, 3.0, NaN], ones(4))
1.0
```
"""
function calc_weighted_std(vals::AbstractArray, weights; dims=:)
    keep, w = _prepare(vals, nothing, weights, dims)
    s1 = _first_moments(vals, nothing, w, keep)
    mx = _ratio.(s1.wx, s1.w)
    s2 = _second_moments(vals, nothing, w, keep, mx, nothing)
    RT = _rtype(vals, nothing, w)
    denom = s1.w .^ 2 .- s2.w2
    variance = ifelse.(denom .> 0, s1.w .* s2.xx ./ denom, NaN)
    return _finish(dims, RT.(sqrt.(max.(variance, 0))))
end

"""
    calc_weighted_pattern_corr(vals1, vals2, weights; dims=:)

Return the centered weighted pattern correlation (Pearson's r) of `vals1` and `vals2` over
`dims`, `Σ w (x - μx)(y - μy) / sqrt(Σ w (x - μx)² Σ w (y - μy)²)`. Only points valid in
both fields are used, also for the means.

# Arguments
- `vals1`, `vals2::AbstractArray{<:Union{Real,Missing}}`: fields of the same size.
- `weights`: weights broadcast-compatible with the fields (an array or a `Real`).

# Keywords
- `dims=:`: dimensions to reduce; `:` returns a scalar.

# Returns
- The correlation in `[-1, 1]`, shaped as for [`calc_weighted_mean`](@ref); `NaN` where a
    field has zero variance or no point is valid.

See also [`calc_correlation_significance`](@ref).

# Examples
```jldoctest
julia> x = [1.0, 2.0, 3.0, 4.0];

julia> calc_weighted_pattern_corr(x, [2.0, 4.0, 6.0, NaN], ones(4))
1.0

julia> calc_weighted_pattern_corr(x, [1.0, 3.0, 2.0, 4.0], [1, 1, 1, 1])
0.8
```
"""
function calc_weighted_pattern_corr(
    vals1::AbstractArray, vals2::AbstractArray, weights; dims=:
)
    keep, w = _prepare(vals1, vals2, weights, dims)
    s1 = _first_moments(vals1, vals2, w, keep)
    mx, my = _ratio.(s1.wx, s1.w), _ratio.(s1.wy, s1.w)
    s2 = _second_moments(vals1, vals2, w, keep, mx, my)
    RT = _rtype(vals1, vals2, w)
    denom = sqrt.(s2.xx .* s2.yy)
    r = ifelse.(denom .> 0, clamp.(s2.xy ./ denom, -1, 1), NaN)
    return _finish(dims, RT.(r))
end

"""
    fisher_z_transform(r)

Return the Fisher z-transform `atanh(r)` of correlation coefficient(s) `r ∈ [-1, 1]`
(elementwise for arrays); `±1` give `±Inf`.

See also [`fisher_z_inv`](@ref).

# Examples
```jldoctest
julia> fisher_z_transform(0.5) ≈ 0.5 * log(3)
true
```
"""
fisher_z_transform(r::Number) = atanh(r)
fisher_z_transform(r::AbstractArray) = atanh.(r)

"""
    fisher_z_inv(z)

Return the correlation coefficient(s) `tanh(z)` for Fisher z-value(s) `z` (elementwise for
arrays), the inverse of [`fisher_z_transform`](@ref).

# Examples
```jldoctest
julia> fisher_z_inv(fisher_z_transform(0.3)) ≈ 0.3
true
```
"""
fisher_z_inv(z::Number) = tanh(z)
fisher_z_inv(z::AbstractArray) = tanh.(z)

"""
    calc_correlation_significance(r, n_eff) -> (p_value, z_score)

Return the two-sided p-value and z-score of correlation coefficient(s) `r` against zero
correlation, from the Fisher z-transform with standard error `1 / sqrt(n_eff - 3)`
[Fisher, 1921]: `z_score = atanh(r) sqrt(n_eff - 3)`, `p_value = erfc(|z_score| / √2)`.

# Arguments
- `r`: correlation coefficient (`Real`) or array of them.
- `n_eff::Real`: effective number of independent samples, `> 3`; for spatial fields
    usually much smaller than the number of grid cells.

# Returns
- `Tuple`: `(p_value, z_score)`, scalars or arrays like `r`.

# Throws
- `ArgumentError`: if `n_eff ≤ 3`.

# References
- Fisher, R. A. (1921), On the "probable error" of a coefficient of correlation deduced
    from a small sample, *Metron*, 1, 3–32. No DOI.

# Examples
```jldoctest
julia> p, z = calc_correlation_significance(0.5, 30);

julia> round(z; digits=3), round(p; digits=5)
(2.854, 0.00431)
```
"""
function calc_correlation_significance(r, n_eff::Real)
    n_eff > 3 || throw(ArgumentError("n_eff must be greater than 3, got $n_eff"))
    z_score = fisher_z_transform(r) .* sqrt(n_eff - 3)
    p_value = erfc.(abs.(z_score) ./ sqrt(2))
    return p_value, z_score
end

"""
    calc_latitude_weights(lat)

Return area weights `cos(lat)` for latitudes `lat` in degrees, normalized to sum 1 and
clipped at 0 beyond the poles. They are proportional to the cell areas of a regular grid
with cell-centered latitudes; for other grids (e.g. half-size polar cells), use the cell
areas instead.

# Arguments
- `lat::AbstractVector{<:Real}`: latitudes in degrees.

# Returns
- `Vector`: normalized weights in the floating-point type of `lat`.

# Examples
```jldoctest
julia> calc_latitude_weights([-60, 0, 60])
3-element Vector{Float64}:
 0.25
 0.5
 0.25
```
"""
function calc_latitude_weights(lat::AbstractVector{<:Real})
    T = float(eltype(lat))
    w = max.(zero(T), cosd.(T.(lat)))
    total = sum(w)
    return total > 0 ? w ./ total : w
end

# --- Internal helpers ---

_isvalid(x) = !(ismissing(x) || isnan(x))
_ratio(a, b) = b > 0 ? a / b : oftype(a / one(b), NaN)

# Result type: floating type of the inputs (missing excluded).
function _rtype(x, y, w)
    T = float(promote_type(nonmissingtype(eltype(x)), eltype(w)))
    return y === nothing ? T : float(promote_type(T, nonmissingtype(eltype(y))))
end

# Accumulator type: at least Float64.
_acctype(x, y, w) = promote_type(_rtype(x, y, w), Float64)

# Check sizes and dims; return the kept dimensions (true = not reduced) and the weights as
# an array.
function _prepare(x, y, weights, dims)
    Base.require_one_based_indexing(x)
    if y !== nothing && size(y) != size(x)
        throw(DimensionMismatch("fields have sizes $(size(x)) and $(size(y))"))
    end
    w = weights isa AbstractArray ? weights : fill(weights)
    Base.require_one_based_indexing(w)
    N = ndims(x)
    if ndims(w) > N || any(d -> size(w, d) ∉ (1, size(x, d)), 1:ndims(w))
        throw(
            DimensionMismatch("weights of size $(size(w)) do not broadcast to $(size(x))")
        )
    end
    if dims === (:)
        keep = ntuple(_ -> false, N)
    else
        ds = dims isa Integer ? (dims,) : Tuple(dims)
        all(d -> 1 <= d <= N, ds) ||
            throw(ArgumentError("invalid dims $dims for $N-d data"))
        keep = ntuple(d -> !(d in ds), N)
    end
    return keep, w
end

_finish(dims, a) = dims === (:) ? first(a) : a

# Index of the reduced output (kept dimensions) and of the broadcast weights.
@inline _project(I::CartesianIndex{N}, keep) where {N} =
    CartesianIndex(ntuple(d -> keep[d] ? I[d] : 1, Val(N)))
@inline _wproject(I::CartesianIndex, wfull::NTuple{M,Bool}) where {M} =
    CartesianIndex(ntuple(d -> wfull[d] ? I[d] : 1, Val(M)))

# First pass: Σw, Σwx, and Σwy over points valid in x (and y) with w > 0.
function _first_moments(x, y, w, keep)
    AT = _acctype(x, y, w)
    osz = ntuple(d -> keep[d] ? size(x, d) : 1, ndims(x))
    sw, sx = zeros(AT, osz), zeros(AT, osz)
    sy = y === nothing ? nothing : zeros(AT, osz)
    wfull = ntuple(d -> size(w, d) != 1, ndims(w))
    @inbounds for I in CartesianIndices(x)
        xi, wi = x[I], w[_wproject(I, wfull)]
        (_isvalid(xi) && wi > 0) || continue
        y === nothing || _isvalid(y[I]) || continue
        J = _project(I, keep)
        sw[J] += wi
        sx[J] += wi * AT(xi)
        y === nothing || (sy[J] += wi * AT(y[I]))
    end
    return (; w=sw, wx=sx, wy=sy)
end

# Second pass: Σw(x - μx)², Σw(y - μy)², Σw(x - μx)(y - μy), and Σw².
function _second_moments(x, y, w, keep, mx, my)
    AT = eltype(mx)
    sxx, sw2 = zeros(AT, size(mx)), zeros(AT, size(mx))
    syy = y === nothing ? nothing : zeros(AT, size(mx))
    sxy = y === nothing ? nothing : zeros(AT, size(mx))
    wfull = ntuple(d -> size(w, d) != 1, ndims(w))
    @inbounds for I in CartesianIndices(x)
        xi, wi = x[I], w[_wproject(I, wfull)]
        (_isvalid(xi) && wi > 0) || continue
        y === nothing || _isvalid(y[I]) || continue
        J = _project(I, keep)
        dx = AT(xi) - mx[J]
        sxx[J] += wi * dx^2
        sw2[J] += AT(wi)^2
        if y !== nothing
            dy = AT(y[I]) - my[J]
            syy[J] += wi * dy^2
            sxy[J] += wi * dx * dy
        end
    end
    return (; xx=sxx, yy=syy, xy=sxy, w2=sw2)
end

end # module StatisticUtils
