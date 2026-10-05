module StatisticUtils

export calc_weighted_pattern_corr, calc_weighted_mean, calc_weighted_std, fisher_z_transform
export fisher_z_inv, calc_correlation_significance, calc_latitude_weights

using Distributions, Statistics

"""
    calc_weighted_pattern_corr(vals1::AbstractArray{T}, vals2::AbstractArray{S}, weights::AbstractArray{W}; dims=:) where {T, S, W}

Calculates the weighted pattern correlation coefficient (Pearson R) between two fields.
Highly optimized for GCM datasets, supporting any rank (e.g., 5-D/6-D) and any floating-point precision.

# Arguments
- `vals1`, `vals2`: Data tensors of identical shape. `NaN` values in either field are masked.
- `weights`: Weighting factors (e.g., grid-cell area). Must be broadcast-compatible with inputs.
- `dims`: Dimensions to reduce. Use `:` (default) for a global scalar or a `Tuple` (e.g., `dims=(1,2)`) for spatial correlation at each time/level.

# Returns
- A scalar or array in promoted precision (`Float32`, `Float64`, `BigFloat`, etc.).
- Returns `NaN` if variance is zero or data is insufficient.

# Method
Uses the weighted Pearson correlation formula:
R = Σ [w * (x - μx) * (y - μy)] / sqrt(Σ [w * (x - μx)²] * Σ [w * (y - μy)²])
"""
function calc_weighted_pattern_corr(vals1::AbstractArray{T}, vals2::AbstractArray{S}, weights::AbstractArray{W}; dims=:) where {T, S, W}

    # Determine promoted result type (e.g., Float32 + Float64 -> Float64)
    RT = promote_type(T, S, W)

    # Calculate weighted means using the generic mean function
    mx = calc_weighted_mean(vals1, weights; dims=dims)
    my = calc_weighted_mean(vals2, weights; dims=dims)

    # Create a unified mask (Valid = No NaN in either field AND weight > 0)
    valid_mask = .!isnan.(vals1) .& .!isnan.(vals2) .& (weights .> zero(W))

    # Fused Broadcasting for Numerator (Covariance) and Denominator (Variances)
    # This prevents the allocation of large intermediate dx and dy arrays
    num_comp = broadcast((v1, v2, w, m, ax, ay) -> m ? RT(w) * (RT(v1) - RT(ax)) * (RT(v2) - RT(ay)) : zero(RT),
                         vals1, vals2, weights, valid_mask, mx, my)

    var_x_comp = broadcast((v1, w, m, ax) -> m ? RT(w) * (RT(v1) - RT(ax))^2 : zero(RT),
                           vals1, weights, valid_mask, mx)

    var_y_comp = broadcast((v2, w, m, ay) -> m ? RT(w) * (RT(v2) - RT(ay))^2 : zero(RT),
                           vals2, weights, valid_mask, my)

    # Reduction Sums
    numerator = sum(num_comp, dims=dims)
    var_x     = sum(var_x_comp, dims=dims)
    var_y     = sum(var_y_comp, dims=dims)

    # Final correlation with zero-variance protection
    denom = @. sqrt(var_x * var_y)
    res = ifelse.(denom .> zero(RT), numerator ./ denom, RT(NaN))

    return dims === (:) ? first(res) : res
end # calc_weighted_pattern_corr

"""
    calc_weighted_mean(vals::AbstractArray{T}, weights::AbstractArray{S}; dims=:) where {T, S}

Calculates the weighted arithmetic mean of a multi-dimensional array, specifically
optimized for planetary climate data (e.g., area-weighted or mass-weighted diagnostics).

# Arguments
- `vals::AbstractArray{T}`: An N-dimensional tensor (e.g., 5-D or 6-D) containing
  the data values. Supports arbitrary precision (`Float16` through `BigFloat`).
  `NaN` values are identified and automatically masked out of the summation.
- `weights::AbstractArray{S}`: Weighting factors that must be broadcast-compatible
  with `vals`. Only elements where `weights > 0` are included.
- `dims`: The dimension(s) over which to perform the reduction.
    - `:` (default): Global reduction to a single scalar.
    - `Int` or `Tuple{Int}`: Reduces specific axes (e.g., `dims=(3,4)` for a
      horizontal mean over lon/lat, preserving vertical and time axes).

# Returns
- A scalar or array in the promoted precision of `T` and `S`. If a reduction
  subset contains no valid data or zero total weight, the result is `NaN`
  in the appropriate precision.

# Numerical Features
- **Precision Preservation**: Uses `promote_type` to ensure intermediate
  calculations occur at the highest precision of the inputs.
- **Memory Efficiency**: Uses broadcast fusion to calculate weighted sums
  in-place, avoiding large intermediate array allocations.
"""
function calc_weighted_mean(vals::AbstractArray{T}, weights::AbstractArray{S}; dims=:) where {T, S}

    # Determine promoted result type (e.g., Float32 + Float64 -> Float64)
    RT = promote_type(T, S)

    # Create a mask for valid data using the weight's zero type
    valid_mask = .!isnan.(vals) .& (weights .> zero(S))

    # Compute weighted values and sum weights using fused broadcasting
    # RT(v) ensures promotion to the higher precision before multiplication
    weighted_vals = broadcast((v, w, m) -> m ? RT(v) * RT(w) : zero(RT), vals, weights, valid_mask)
    valid_weights = broadcast((w, m) -> m ? RT(w) : zero(RT), weights, valid_mask)

    # Perform the reduction
    sum_wv = sum(weighted_vals, dims=dims)
    sum_w  = sum(valid_weights, dims=dims)

    # Final division with NaN protection using type-appropriate NaN
    res = ifelse.(sum_w .> zero(RT), sum_wv ./ sum_w, RT(NaN))

    return dims === (:) ? first(res) : res
end # calc_weighted_mean

"""
    calc_weighted_std(vals::AbstractArray{T}, weights::AbstractArray{S}; dims=:) where {T, S}

Calculates the weighted standard deviation using the **Reliability Weights** correction
(unbiased estimator for weighted data). This is essential for  calculating spatial or
temporal variability in GCM outputs where grid cells have unequal areas or masses.

# Arguments
- `vals::AbstractArray{T}`: Data tensor of arbitrary rank and precision.
  `NaN` values are ignored.
- `weights::AbstractArray{S}`: Weights (e.g., area or pressure) broadcast-compatible
  with `vals`. Non-positive weights are ignored.
- `dims`: Dimensions to collapse. Default is `:` (global reduction).

# Returns
- The weighted standard deviation in the promoted precision of the inputs.
  Returns `NaN` if the effective sample size is insufficient (denominator <= 0).

# Mathematical Logic
```math
\\sigma_w = \\sqrt{ \\frac{V_1 \\sum_{i=1}^N w_i(x_i - \\mu_w)^2}{V_1^2 - V_2} }
"""
function calc_weighted_std(vals::AbstractArray{T}, weights::AbstractArray{S}; dims=:) where {T, S}

    # Determine promoted result type (e.g., Float32 + Float64 -> Float64)
    RT = promote_type(T, S)

    # Calculate weighted mean (automatically inherits the promoted type RT)
    mx = calc_weighted_mean(vals, weights; dims=dims)

    # Define valid mask
    valid_mask = .!isnan.(vals) .& (weights .> zero(S))

    # Fused broadcasting for variance components (V1 and V2)
    # This pass calculates the squared differences in-place
    weighted_sq_diff = broadcast((v, w, m, avg) -> m ? RT(w) * (RT(v) - RT(avg))^2 : zero(RT),
                                 vals, weights, valid_mask, mx)

    w_sum   = broadcast((w, m) -> m ? RT(w) : zero(RT), weights, valid_mask)
    w2_sum  = broadcast((w, m) -> m ? RT(w)^2 : zero(RT), weights, valid_mask)

    # Reduction
    sum_wsd = sum(weighted_sq_diff, dims=dims)
    V1      = sum(w_sum, dims=dims)
    V2      = sum(w2_sum, dims=dims)

    # Reliability weight correction: (V1 * sum_wsd) / (V1^2 - V2)
    denom = @. V1^2 - V2
    variance = ifelse.(denom .> zero(RT), (sum_wsd .* V1) ./ denom, RT(NaN))

    # Final sqrt with floor at zero to prevent imaginary numbers from roundoff
    res = sqrt.(max.(zero(RT), variance))

    return dims === (:) ? first(res) : res
end # calc_weighted_std

"""
    fisher_z_transform(r::AbstractArray{T}) where {T}
    fisher_z_transform(r::Number)

Converts Pearson correlation coefficients (r) to Fisher Z-values using:
Z = 0.5 * ln((1 + r) / (1 - r)) or Z = atanh(r).

# Arguments
- `r`: Correlation coefficient(s) in the range (-1, 1).

# Returns
- Z-transformed values in the same precision as the input.
  Values of exactly 1.0 or -1.0 return Inf or -Inf.
"""
function fisher_z_transform(r::AbstractArray{T}) where {T}
    RT = eltype(r)
    # Clamp r slightly to avoid Log(0) if r is exactly 1 or -1
    return atanh.(r)
end
fisher_z_transform(r::Number) = atanh(r)

"""
    fisher_z_inv(z::AbstractArray{T}) where {T}
    fisher_z_inv(z::Number)

Converts Fisher Z-values back to Pearson correlation coefficients (r) using:
r = (exp(2Z) - 1) / (exp(2Z) + 1) or r = tanh(Z).
"""
function fisher_z_inv(z::AbstractArray{T}) where {T}
    return tanh.(z)
end
fisher_z_inv(z::Number) = tanh(z)

"""
    calc_correlation_significance(r, n_eff)

Calculates the p-value and Z-score for a given pattern correlation `r`
using the Fisher Z-transformation and the Effective Degrees of Freedom `n_eff`.

# Arguments
- `r`: Pearson correlation coefficient (scalar or array).
- `n_eff`: The Effective Degrees of Freedom (EDOF). For spatial fields,
  this is typically much smaller than the total number of grid cells.

# Returns
- `p_value`: The two-tailed probability that the correlation is due to chance.
- `z_score`: The number of standard deviations the transformed Z is from zero.

# Mathematical Method
1. Transform to Z: ``Z = atanh(r)``
2. Standard Error: ``SE = 1 / sqrt(n_eff - 3)``
3. Z-score: ``Z_{score} = Z / SE``
"""
function calc_correlation_significance(r::AbstractArray{T}, n_eff::Number) where {T}

    # Fisher Z-transform
    z_vals = atanh.(r)

    # Calculate Standard Error based on EDOF
    # Note: n_eff must be > 3 for a valid SE
    se = 1.0 / sqrt(n_eff - 3.0)

    # Calculate Z-score (how many standard errors away from 0)
    z_scores = z_vals ./ se

    # Calculate two-tailed p-value using the Standard Normal Distribution
    dist = Normal(0, 1)
    p_values = 2 .* (1 .- cdf.(Ref(dist), abs.(z_scores)))

    return p_values, z_scores
end # calc_correlation_significance

"""
    calc_latitude_weights(lat::AbstractVector{T}) where {T}

Calculates spatial weights based on the cosine of latitude for a regular lat/lon grid.
Weights represent the relative surface area of grid cells at different latitudes.

# Arguments
- `lat::AbstractVector{T}`: Vector of latitudes in degrees (-90 to 90).
  Supports any precision (Float32, Float64, BigFloat).

# Returns
- A vector of normalized weights (sum = 1.0) of the same length and type as `lat`.

# Method
1. Converts latitude from degrees to radians.
2. Calculates the cosine of the latitude: weight = cos(lat * pi / 180).
3. Ensures all weights are non-negative and normalizes the sum to 1.0.
"""
function calc_latitude_weights(lat::AbstractVector{T}) where {T}

    # Convert degrees to radians and calculate cosine
    # Use T(pi) to match the precision of the input array
    weights = cos.(lat .* (T(π) / T(180)))

    # Ensure weights are non-negative (clipping at poles)
    weights .= max.(zero(T), weights)

    # Normalize weights so they sum to 1.0
    total_w = sum(weights)
    if total_w > zero(T)
        weights ./= total_w
    end

    return weights
end
end # module StatisticUtils
