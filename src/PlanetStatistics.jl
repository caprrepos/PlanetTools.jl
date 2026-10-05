module PlanetStatistics

export calc_planet_means!, calc_annual_mean_ts!, calc_seasonal_mean_ts!, calc_multiyear_stats!
export calc_anomalies!, calc_linear_trends!, calc_linear_trends_sig!, generate_significance_mask!

using Base.Threads
using Dates
using Distributions
using NaNStatistics
using Statistics


# --- Private Helper Functions ---

"""
    _nansum(v; dims)

Calculates the sum of array `v` along dimensions `dims`, strictly treating `NaN`
values as zero. This ensures numerical stability in planetary datasets where
missing values (e.g., polar night gaps, masked topography, or sensor dropouts)
are expected.

### Mathematical Definition
For an N-dimensional slice \$\\mathbf{x}\$:
\$\$ \\sum_{i} x_i \\quad \\forall x_i \\in \\mathbb{R} \\setminus \\{\\text{NaN}\\} \$\$

Returns the additive identity (zero) of the input's element type if an entire
slice is composed of `NaN`s.
"""
function _nansum(v; dims)
    return sum(ifelse.(isnan.(v), zero(eltype(v)), v); dims=dims)
end

"""
    calc_planet_means!(data, zm_out, gm_out, weights_zm, weights_gm, axyp, t_ranges;
                       run_names=nothing, threshold=100.0)

Calculates the zonal and area-weighted global means for planetary fields across multiple
simulations, while tracking the fraction of valid data points.

### Arguments
- `data`: Input array with shape `(Sim, Var, Lon, Lat, [Lev], Time)`.
- `zm_out`: Output array for zonal means: `(Sim, Var, 1, Lat, [Lev], Time)`.
- `gm_out`: Output array for global means: `(Sim, Var, 1, 1, [Lev], Time)`.
- `weights_zm`: Output array for zonal availability fraction: `(Sim, Var, 1, Lat, [Lev], Time)`.
- `weights_gm`: Output array for global availability fraction: `(Sim, Var, 1, 1, [Lev], Time)`.
- `axyp`: 2D array of grid cell areas `(Lon, Lat)` used for global area-weighting.
- `t_ranges`: Vector of integers specifying the time length of each simulation.
- `run_names`: (Optional) Vector of strings for logging.
- `threshold`: (Optional) Minimum percentage (0-100) of non-NaN grid points required.

### Logic
1. **Zonal Mean**: Averages across Lon (Dim 3). Weights represent the fraction of valid longitudes.
2. **Global Mean**: Area-weighted average across Lon/Lat (Dims 3 & 4). Weights represent
   the fraction of the horizontal plane containing valid data.
3. **Thresholding**: If availability < `threshold`, output data is `NaN` and weight is `0.0`.
   Otherwise, weight is `valid_count / total_possible_count`.

# Returns
- `nothing`
"""
function calc_planet_means!(data, zm_out, gm_out, weights_zm, weights_gm, axyp, t_ranges;
                            run_names=nothing, threshold=100.0)
    nd = ndims(data)
    LON_DIM = 3
    LAT_DIM = 4

    # Pre-calculate area weights for broadcasting (1, 1, Lon, Lat, 1, 1)
    weight_dims = ntuple(i -> (i == 3 || i == 4) ? size(axyp, i-2) : 1, nd)
    axyp_nd = reshape(Float32.(axyp), weight_dims)

    # Dimensional counts for fraction calculation
    n_lon = size(data, LON_DIM)
    n_horiz = size(data, LON_DIM) * size(data, LAT_DIM)

    # Parallel loop over simulation indices
    Threads.@threads for ns in eachindex(t_ranges)
        name = isnothing(run_names) ? "Sim $ns" : run_names[ns]
        println("Thread $(Threads.threadid()) - Zonal and global means for run: ", name)

        t_limit = t_ranges[ns]
        t_range = 1:t_limit

        # Slice simulation-specific data
        colons = ntuple(_ -> :, nd - 2)
        data_view = @view data[ns:ns, colons..., t_range]

        # Identify valid (non-NaN) data presence
        mask = .!isnan.(data_view)

        # --- Zonal Mean and Weights ---
        # Calculate fraction of non-missing longitudes
        valid_count_zm = sum(mask, dims=LON_DIM)
        fraction_zm = Float32.(valid_count_zm ./ n_lon)
        presence_zm = fraction_zm .* 100.0f0

        zm_raw = nanmean(data_view, dims=LON_DIM)

        # Apply threshold to data and weights
        zm_out[ns:ns, colons..., t_range] .= ifelse.(presence_zm .< threshold,
                                                    Float32(NaN), zm_raw)
        weights_zm[ns:ns, colons..., t_range] .= ifelse.(presence_zm .< threshold,
                                                        0.0f0, fraction_zm)

        # --- Weighted Global Mean and Weights ---
        # Calculate fraction of non-missing horizontal grid points
        valid_count_gm = sum(mask, dims=(LON_DIM, LAT_DIM))
        fraction_gm = Float32.(valid_count_gm ./ n_horiz)
        presence_gm = fraction_gm .* 100.0f0

        # Numerator: sum(data * area) | Denominator: sum(area) where data exists
        num = _nansum(data_view .* axyp_nd, dims=(LON_DIM, LAT_DIM))
        den = sum(mask .* axyp_nd, dims=(LON_DIM, LAT_DIM))
        gm_raw = num ./ den

        # Apply threshold to data and weights
        gm_out[ns:ns, colons..., t_range] .= ifelse.(presence_gm .< threshold,
                                                    Float32(NaN), gm_raw)
        weights_gm[ns:ns, colons..., t_range] .= ifelse.(presence_gm .< threshold,
                                                        0.0f0, fraction_gm)
    end

    return nothing
end # function calc_planet_means!

"""
    calc_annual_mean_ts!(data, zm, gm,
                         ann_out, zm_ann_out, gm_ann_out,
                         ann_w_out, zm_ann_w_out, gm_ann_w_out,
                         dates, t_ranges, fullYearDates;
                         data_w=nothing, zm_w=nothing, gm_w=nothing,
                         month_weights=nothing, run_names=nothing,
                         threshold=100.0)

Calculates a month-weighted annual mean time series for grid-point, zonal-mean, and
global-mean fields, integrating optional spatial availability weights and enforcing
data-sufficiency thresholds.

# Arguments

## Primary Input Data Arrays
- `data`: Raw field data. Shape: `(Sim, Var, Lon, Lat, [Lev], Time)`.
- `zm`: Pre-calculated zonal means. Shape: `(Sim, Var, 1, Lat, [Lev], Time)`.
- `gm`: Pre-calculated global means. Shape: `(Sim, Var, 1, 1, [Lev], Time)`.

## Optional Input Weight Arrays
- `data_w`: Optional spatial weights for `data`. Shape must match `data`.
- `zm_w`: Optional availability weights for `zm`. Shape must match `zm`.
- `gm_w`: Optional availability weights for `gm`. Shape must match `gm`.
  *Note: These typically represent the fraction of valid data within a spatial reduction
  (e.g., from `calc_planet_means!`). If `nothing`, defaults to 1.0 for non-NaN data.*

## Output Mean Arrays (Modified In-Place)
- `ann_out`: Annual grid-point means. Shape: `(Sim, Var, Lon, Lat, [Lev], Year)`.
- `zm_ann_out`: Annual zonal means. Shape: `(Sim, Var, 1, Lat, [Lev], Year)`.
- `gm_ann_out`: Annual global means. Shape: `(Sim, Var, 1, 1, [Lev], Year)`.

## Output Weight Arrays (Modified In-Place)
- `ann_w_out`: Annual grid-point weights. Shape: `(Sim, Var, Lon, Lat, [Lev], Year)`.
- `zm_ann_w_out`: Annual zonal weights. Shape: `(Sim, Var, 1, Lat, [Lev], Year)`.
- `gm_ann_w_out`: Annual global weights. Shape: `(Sim, Var, 1, 1, [Lev], Year)`.
  *Note: These store the final annual availability fraction [0.0, 1.0] per year.*

## Metadata & Control
- `dates`: `Vector{Date}` representing the master timeline of the input arrays.
- `t_ranges`: `Vector{Int}` containing the number of monthly steps per simulation.
- `fullYearDates`: `Vector{StepRange{Date}}` defining target years for each simulation.
- `month_weights`: (Optional) `Vector` of temporal weights (e.g., days in month).
- `run_names`: (Optional) `Vector{String}` for logging simulation progress.
- `threshold`: (Optional) Minimum percentage (0-100) of valid months required per year.

# Logic
1. **Combined Weighting**: For each month, the total weight is calculated as:
   `Combined_Weight = Input_Weight * Month_Weight`.
2. **Availability Calculation**:
   - `Fraction = (Count of months where Input_Weight > 0 and Data is not NaN) / 12`.
   - If `(Fraction * 100) < threshold`, the output weight is `0.0` and the mean is `NaN`.
   - Otherwise, the output weight is the calculated fraction.
3. **Normalization**: The annual average is computed by dividing the sum of weighted
   monthly values by the sum of the `Combined_Weights` for valid months.

# Returns
- `nothing`
"""
function calc_annual_mean_ts!(data, zm, gm,
                              ann_out, zm_ann_out, gm_ann_out,
                              ann_w_out, zm_ann_w_out, gm_ann_w_out,
                              dates, t_ranges, fullYearDates;
                              data_w=nothing, zm_w=nothing, gm_w=nothing,
                              month_weights=nothing, run_names=nothing,
                              threshold=100.0)

    nd_data = ndims(data)
    nd_zm   = ndims(zm)
    nd_gm   = ndims(gm)
    MONTHS_IN_YEAR = 12

    # Master temporal weights (e.g., days in month)
    m_weights = isnothing(month_weights) ? ones(Float32, length(dates)) : month_weights

    Threads.@threads for ns in eachindex(t_ranges)
        name = isnothing(run_names) ? "Sim $ns" : run_names[ns]
        println("Thread $(Threads.threadid()) - Annual means for run: ", name)

        sim_years_range = fullYearDates[ns]
        if isempty(sim_years_range) continue end

        target_start = first(sim_years_range)
        nd_start = findfirst(d -> year(d) == year(target_start) &&
                                 month(d) == month(target_start), dates)

        if isnothing(nd_start)
            @warn "Start date for $name not found in master dates vector."
            continue
        end

        # Reshape month_weights for broadcasting (1, 1, ..., Time)
        mw_ij = reshape(m_weights, ntuple(_ -> 1, nd_data - 2)..., :)
        mw_zm = reshape(m_weights, ntuple(_ -> 1, nd_zm - 2)..., :)
        mw_gm = reshape(m_weights, ntuple(_ -> 1, nd_gm - 2)..., :)

        colons_data = ntuple(_ -> :, nd_data - 2)
        colons_zm   = ntuple(_ -> :, nd_zm - 2)
        colons_gm   = ntuple(_ -> :, nd_gm - 2)

        curr_nd = nd_start

        for (ny, _) in enumerate(sim_years_range)
            t_slice = curr_nd:(curr_nd + MONTHS_IN_YEAR - 1)
            if last(t_slice) > t_ranges[ns] break end

            # Internal engine for combined weighting and thresholding
            function weighted_stats_combined(d_v, in_w_v, m_w_v, dim, thresh)
                # Default input weight to 1.0/0.0 based on NaN if not provided
                eff_in_w = isnothing(in_w_v) ? ifelse.(isnan.(d_v), 0.0f0, 1.0f0) : in_w_v

                # Combined weight = spatial availability * temporal length
                combined_w = eff_in_w .* m_w_v
                mask = .!isnan.(d_v) .& (eff_in_w .> 0.0f0)

                valid_count = sum(mask, dims=dim)
                fraction = Float32.(valid_count ./ MONTHS_IN_YEAR)
                presence_pct = fraction .* 100.0f0

                num = nansum(d_v .* combined_w, dims=dim)
                den = sum(combined_w, dims=dim)
                res_mean = num ./ den

                is_invalid = presence_pct .< thresh
                final_mean = ifelse.(is_invalid, eltype(res_mean)(NaN), res_mean)
                final_weight = ifelse.(is_invalid, 0.0f0, fraction)

                return final_mean, final_weight
            end

            # --- Grid-Point Processing ---
            data_v = @view data[ns, colons_data..., t_slice]
            in_w_v = isnothing(data_w) ? nothing : @view data_w[ns, colons_data..., t_slice]
            m_ij, w_ij = weighted_stats_combined(data_v, in_w_v, @view(mw_ij[colons_data..., t_slice]),
                nd_data - 1, threshold)
            ann_out[ns, colons_data..., ny] .= dropdims(m_ij, dims=nd_data-1)
            ann_w_out[ns, colons_data..., ny] .= dropdims(w_ij, dims=nd_data-1)

            # --- Zonal Processing ---
            zm_v = @view zm[ns, colons_zm..., t_slice]
            in_zw_v = isnothing(zm_w) ? nothing : @view zm_w[ns, colons_zm..., t_slice]
            m_zm, w_zm = weighted_stats_combined(zm_v, in_zw_v, @view(mw_zm[colons_zm..., t_slice]),
                nd_zm - 1, threshold)
            zm_ann_out[ns, colons_zm..., ny] .= dropdims(m_zm, dims=nd_zm-1)
            zm_ann_w_out[ns, colons_zm..., ny] .= dropdims(w_zm, dims=nd_zm-1)

            # --- Global Processing ---
            gm_v = @view gm[ns, colons_gm..., t_slice]
            in_gw_v = isnothing(gm_w) ? nothing : @view gm_w[ns, colons_gm..., t_slice]
            m_gm, w_gm = weighted_stats_combined(gm_v, in_gw_v, @view(mw_gm[colons_gm..., t_slice]),
                nd_gm - 1, threshold)
            gm_ann_out[ns, colons_gm..., ny] .= dropdims(m_gm, dims=nd_gm-1)
            gm_ann_w_out[ns, colons_gm..., ny] .= dropdims(w_gm, dims=nd_gm-1)

            curr_nd += MONTHS_IN_YEAR
        end
    end

    return nothing
end # function calc_annual_mean_ts!

"""
    calc_weighted_seasonal_mean(data_slice, in_w_slice, m_w_slice, threshold, length_of_season)

Helper engine that performs a weighted temporal reduction on a multi-dimensional seasonal
data slice, accounting for both spatial and temporal weights.

# Arguments
- `data_slice`: A view of the data for a specific season.
  Shape: `(Var, [Spatial...], Time)`.
- `in_w_slice`: Optional spatial/fractional weights (e.g., from zonal/global mean steps).
  Shape: Must match `data_slice`. If `nothing`, defaults to `1.0` for non-NaN data.
- `m_w_slice`: 1D vector of temporal weights (e.g., days per month) for the season.
  Length: Must match the `Time` dimension of `data_slice`.
- `threshold`: Minimum percentage (0.0 to 100.0) of valid months required.
- `length_of_season`: Expected number of months in the season (e.g., 3 for DJF).

# Mathematical Logic
1. Total Weighting: Combined weight W = input_weight * month_weight.
2. Mean: avg = sum(x * W) / sum(W) for all non-NaN x where W > 0.
3. Availability: Fraction f = (count of valid months) / length_of_season.
4. Sufficiency: If (f * 100) < threshold, Mean is NaN and Weight is 0.0.

# Returns
- `(mean_res, weight_res)`: A tuple of arrays with rank \$N-1\$ (Time dimension dropped).
"""
function _calc_weighted_seasonal_mean(data_slice, in_w_slice, m_w_slice, threshold,
        length_of_season)
    t_dim = ndims(data_slice)
    el_type = eltype(data_slice)

    # Combined weighting logic
    eff_in_w = isnothing(in_w_slice) ? ifelse.(isnan.(data_slice), 0.0f0, 1.0f0) : in_w_slice
    m_w_b = reshape(m_w_slice, ntuple(_ -> 1, t_dim - 1)..., :)
    combined_w = eff_in_w .* m_w_b

    mask = .!isnan.(data_slice) .& (eff_in_w .> 0.0f0)
    v_count = sum(mask, dims=t_dim)

    # Weighted calculation
    num = _nansum(data_slice .* combined_w, dims=t_dim)
    den = sum(combined_w, dims=t_dim)
    res = num ./ den

    # Availability fraction [0.0 - 1.0]
    fraction = Float32.(v_count ./ length_of_season)
    is_invalid = (fraction .* 100.0f0) .< threshold

    res_m = ifelse.(is_invalid, el_type(NaN), res)
    res_w = ifelse.(is_invalid, 0.0f0, fraction)

    return dropdims(res_m, dims=t_dim), dropdims(res_w, dims=t_dim)
end #  calc_weighted_seasonal_mean

"""
    calc_seasonal_mean_ts!(data, zm, gm, out_ij, out_zm, out_gm,
                           w_out_ij, w_out_zm, w_out_gm,
                           indSeasFirst, dates;
                           data_w=nothing, zm_w=nothing, gm_w=nothing,
                           month_weights=nothing, threshold=100.0,
                           length_of_season=3, run_names=nothing)

Calculates a robust, month-weighted seasonal mean time series for grid-point maps,
zonal-mean sections, and global-mean profiles across multiple simulations.

# Arguments

## Input Data Arrays (Rank Agnostic)
- `data`: Input grid-point data. Shape: `(Sim, Var, Lon, Lat, [Lev], Time)`.
- `zm`: Input zonal-mean data. Shape: `(Sim, Var, 1, Lat, [Lev], Time)`.
- `gm`: Input global-mean data. Shape: `(Sim, Var, 1, 1, [Lev], Time)`.

## Optional Input Weight Arrays
- `data_w`, `zm_w`, `gm_w`: Optional spatial/availability weights.
  Shapes must match their corresponding data arrays exactly.
  Used to weight the months based on spatial data coverage.

## Output Mean Arrays (Modified In-Place)
- `out_ij`: Seasonal grid means. Shape: `(Sim, Var, Season, Lon, Lat, [Lev], Year)`.
- `out_zm`: Seasonal zonal means. Shape: `(Sim, Var, Season, 1, Lat, [Lev], Year)`.
- `out_gm`: Seasonal global means. Shape: `(Sim, Var, Season, 1, 1, [Lev], Year)`.

## Output Weight Arrays (Modified In-Place)
- `w_out_ij`, `w_out_zm`, `w_out_gm`: Final seasonal availability weights.
  Shapes match the Output Mean arrays. Values represent the fraction of valid
  months [0.0, 1.0] used in the seasonal average.

## Metadata & Control
- `indSeasFirst`: `(Sim, Season, Year)` index mapping to the `dates` vector.
- `dates`: Chronological `Vector{Date}` corresponding to the `Time` dimension.
- `month_weights`: (Optional) 1D `Vector` of temporal weights (e.g., days/month).
- `threshold`: (Optional) Min % of valid months required per season. Default 100.0.
- `length_of_season`: (Optional) Months per season. Default 3.
- `run_names`: (Optional) `Vector{String}` for thread-safe logging.

# Returns
- `nothing`
"""
function calc_seasonal_mean_ts!(data, zm, gm, out_ij, out_zm, out_gm,
                                w_out_ij, w_out_zm, w_out_gm,
                                indSeasFirst, dates;
                                data_w=nothing, zm_w=nothing, gm_w=nothing,
                                month_weights=nothing, threshold=100.0,
                                length_of_season=3, run_names=nothing)

    # Identify spatial geometry based on the rank of input tensors
    sp_cols_data = ntuple(_ -> :, ndims(data) - 3)
    sp_cols_zm   = ntuple(_ -> :, ndims(zm) - 3)
    sp_cols_gm   = ntuple(_ -> :, ndims(gm) - 3)

    w_vec = isnothing(month_weights) ? ones(Float32, length(dates)) : month_weights

    Threads.@threads for ns in axes(data, 1)
        name = isnothing(run_names) ? "Sim $ns" : run_names[ns]
        println("Thread $(Threads.threadid()) - Seasonal means for run: ", name)

        for ne in axes(indSeasFirst, 2)     # Seasons
            for ny in axes(indSeasFirst, 3) # Years
                t_start = indSeasFirst[ns, ne, ny]
                if t_start == 0 break end

                t_slice = t_start:(t_start + length_of_season - 1)
                m_w_slice = @view w_vec[t_slice]

                # --- Grid-Point Mean & Weight ---
                d_v = @view data[ns, :, sp_cols_data..., t_slice]
                iw_v = isnothing(data_w) ? nothing : @view data_w[ns, :, sp_cols_data..., t_slice]
                res, w_res = _calc_weighted_seasonal_mean(d_v, iw_v, m_w_slice, threshold,
                    length_of_season)

                out_ij[ns, :, ne, sp_cols_data..., ny] .= res
                w_out_ij[ns, :, ne, sp_cols_data..., ny] .= w_res

                # --- Zonal Mean & Weight ---
                z_v = @view zm[ns, :, sp_cols_zm..., t_slice]
                izw_v = isnothing(zm_w) ? nothing : @view zm_w[ns, :, sp_cols_zm..., t_slice]
                z_res, zw_res = _calc_weighted_seasonal_mean(z_v, izw_v, m_w_slice, threshold,
                    length_of_season)

                out_zm[ns, :, ne, sp_cols_zm..., ny] .= z_res
                w_out_zm[ns, :, ne, sp_cols_zm..., ny] .= zw_res

                # --- Global Mean & Weight ---
                g_v = @view gm[ns, :, sp_cols_gm..., t_slice]
                igw_v = isnothing(gm_w) ? nothing : @view gm_w[ns, :, sp_cols_gm..., t_slice]
                g_res, gw_res = _calc_weighted_seasonal_mean(g_v, igw_v, m_w_slice, threshold,
                    length_of_season)

                out_gm[ns, :, ne, sp_cols_gm..., ny] .= g_res
                w_out_gm[ns, :, ne, sp_cols_gm..., ny] .= gw_res
            end
        end
    end

    return nothing
end # function calc_seasonal_mean_ts!

"""
    _core_statistics!(moa, soa, w_moa, ska, kua, data, time_indices, weights, data_w, threshold)

The internal computational engine for multi-dimensional weighted population statistics.

# Arguments
- `moa::AbstractArray`: **Output** Mean of Analysis buffer. Shape: `(Var, [Spatial...])`.
- `soa::AbstractArray`: **Output** Standard deviation buffer. Shape: `(Var, [Spatial...])`.
- `w_moa::AbstractArray`: **Output** Aggregate weight buffer. Shape: `(Var, [Spatial...])`.
  Stores the fraction of valid time steps [0.0, 1.0].
- `ska`, `kua`: **Output** Skewness and Kurtosis placeholders (currently `nothing`).
- `data::AbstractArray`: **Input** Data view. Last dimension must be the temporal axis.
- `time_indices::AbstractVector{Int}`: Indices to extract from the temporal axis.
- `weights::Union{AbstractVector, Nothing}`: 1D temporal weights (e.g., days per month).
- `data_w::Union{AbstractArray, Nothing}`: **Input** Optional availability weights matching
  `data` shape. If `nothing`, defaults to `1.0` for non-NaN points.
- `threshold::Float64`: Percentage (0-100) of valid data points required to return a result.

# Mathematical Implementation
1. **Total Weighting**: Combined weight \$W_i = data\\_w_i \\times weights_i\$.
2. **Availability**: Fraction \$f = (\\text{count of valid steps}) / (\\text{total window})\$.
3. **Weighted Mean**: \$\\mu = \\sum(x_i W_i) / \\sum(W_i)\$ for valid points.
4. **Weighted Variance**: \$\\sigma^2 = \\sum(W_i (x_i - \\mu)^2) / \\sum(W_i)\$.
5. **Sufficiency**: If \$(f \\times 100) < threshold\$, outputs are set to `NaN` and weight to `0.0`."""
function _core_statistics!(moa, soa, w_moa, ska, kua, data, time_indices, weights, data_w, threshold)
    t_axis = ndims(data)
    v = selectdim(data, t_axis, time_indices)
    window_len = length(time_indices)
    t_dim = ndims(v)

    # 1. Spatial/Availability Weights (data_w)
    # Extract the matching time slice from the input weight tensor
    eff_data_w = isnothing(data_w) ? ifelse.(isnan.(v), 0.0f0, 1.0f0) : selectdim(data_w, t_axis, time_indices)

    # 2. Temporal Weights (e.g. days in month)
    if isnothing(weights)
        w_temp = ones(Float32, ntuple(_ -> 1, t_dim - 1)..., window_len)
    else
        w_temp = reshape(Float32.(weights), ntuple(_ -> 1, t_dim - 1)..., :)
    end

    # 3. Combine weights and mask
    combined_w = eff_data_w .* w_temp
    mask = .!isnan.(v) .& (eff_data_w .> 0.0f0)

    v_count = sum(mask, dims=t_dim)
    sum_w = sum(combined_w, dims=t_dim)

    # Weighted Mean
    mu_w = _nansum(v .* combined_w, dims=t_dim) ./ sum_w

    # Weighted Standard Deviation
    diff = v .- mu_w
    var_w = _nansum(combined_w .* mask .* diff.^2, dims=t_dim) ./ sum_w
    sigma_w = sqrt.(var_w)

    # Calculate Aggregate Availability Weight [0.0 - 1.0]
    fraction = Float32.(v_count ./ window_len)

    # Thresholding
    invalid = dropdims((fraction .* 100.0f0) .< Float64(threshold), dims=t_dim)

    moa .= ifelse.(invalid, eltype(moa)(NaN), dropdims(mu_w, dims=t_dim))
    soa .= ifelse.(invalid, eltype(soa)(NaN), dropdims(sigma_w, dims=t_dim))
    w_moa .= ifelse.(invalid, 0.0f0, dropdims(fraction, dims=t_dim))

    return nothing
end

# -------------------------------------------------------------------------
# 📅 ANNUAL DISPATCH
# -------------------------------------------------------------------------

"""
    calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, indices, Val(:annual); ...)

Calculates multi-year statistics over a continuous block of years.

# Arguments
- `data`: (Sim, Var, [Spatial...], Time).
- `data_w`: (Optional) Weights matching `data` shape.
- `moa`, `soa`, `w_moa`: Output tensors. Shape: (Sim, Var, [Spatial...]).
- `indices`: [Sim, 2] matrix of (start_idx, end_idx).
"""
function calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, indices::AbstractMatrix, ::Val{:annual}; weights=nothing, data_w=nothing, threshold=100.0)
    if !isnothing(data_w) && size(data) != size(data_w)
        throw(DimensionMismatch("Input Data size $(size(data)) and data_w size $(size(data_w)) must match."))
    end

    Threads.@threads for ns in axes(data, 1)
        i_s, i_e = indices[ns, 1], indices[ns, 2]
        if (i_s > 0 && i_e >= i_s)
            _core_statistics!(
                selectdim(moa, 1, ns), selectdim(soa, 1, ns), selectdim(w_moa, 1, ns),
                nothing, nothing,
                selectdim(data, 1, ns), i_s:i_e,
                isnothing(weights) ? nothing : weights[ns, i_s:i_e],
                isnothing(data_w) ? nothing : selectdim(data_w, 1, ns),
                threshold
            )
        end
    end

    return nothing
end

# -------------------------------------------------------------------------
# 📅 MONTHLY DISPATCH
# -------------------------------------------------------------------------

"""
    calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, indices, Val(:monthly); ...)

Calculates climatological monthly statistics (e.g., Average January over 30 years).

# Arguments
- `data`: (Sim, Var, [Spatial...], Time).
- `data_w`: (Optional) Weights matching `data` shape.
- `moa`, `soa`, `w_moa`: Output tensors. Shape: (Sim, Var, [Spatial...], 12).
- `indices`: [Sim, 12] matrix of `Vector{Int}` containing time steps for each month.
"""
function calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, indices::AbstractMatrix, ::Val{:monthly}; weights=nothing, data_w=nothing, threshold=100.0)
    if !isnothing(data_w) && size(data) != size(data_w)
        throw(DimensionMismatch("Input Data size $(size(data)) and data_w size $(size(data_w)) must match."))
    end

    Threads.@threads for ns in axes(data, 1)
        for m in 1:12
            t_vec = indices[ns, m]
            if !isempty(t_vec)
                # Select the specific month in the last dimension of output tensors
                m_view = selectdim(selectdim(moa, 1, ns), ndims(moa)-1, m)
                s_view = selectdim(selectdim(soa, 1, ns), ndims(soa)-1, m)
                w_view = selectdim(selectdim(w_moa, 1, ns), ndims(w_moa)-1, m)

                _core_statistics!(
                    m_view, s_view, w_view, nothing, nothing,
                    selectdim(data, 1, ns), t_vec,
                    isnothing(weights) ? nothing : weights[ns, t_vec],
                    isnothing(data_w) ? nothing : selectdim(data_w, 1, ns),
                    threshold
                )
            end
        end
    end

    return nothing
end

# -------------------------------------------------------------------------
# 📅 SEASONAL DISPATCH
# -------------------------------------------------------------------------

"""
    calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, indices, Val(:seasonal); ...)

Calculates multi-year statistics for seasonal data.

# Arguments
- `data`: (Sim, Var, Season, [Spatial...], Year).
- `data_w`: (Optional) Weights matching `data` shape.
- `moa`, `soa`, `w_moa`: Output tensors. Shape: (Sim, Var, Season, [Spatial...]).
- `indices`: [Sim, Season, 2] matrix of (start_year, end_year).
"""
function calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, indices::AbstractArray{Int, 3}, ::Val{:seasonal}; weights=nothing, data_w=nothing, threshold=100.0)
    if !isnothing(data_w) && size(data) != size(data_w)
        throw(DimensionMismatch("Input Data size $(size(data)) and data_w size $(size(data_w)) must match."))
    end

    Threads.@threads for ns in axes(data, 1)
        for s in axes(indices, 2)
            i_s, i_e = indices[ns, s, 1], indices[ns, s, 2]
            if (i_s > 0 && i_e >= i_s)
                # Seasonal dim is 3. Years dim is last.
                m_view = selectdim(selectdim(moa, 1, ns), 2, s)
                s_view = selectdim(selectdim(soa, 1, ns), 2, s)
                w_view = selectdim(selectdim(w_moa, 1, ns), 2, s)
                d_view = selectdim(selectdim(data, 1, ns), 2, s)
                dw_view = isnothing(data_w) ? nothing : selectdim(selectdim(data_w, 1, ns), 2, s)

                _core_statistics!(
                    m_view, s_view, w_view, nothing, nothing,
                    d_view, i_s:i_e,
                    isnothing(weights) ? nothing : weights[ns, i_s:i_e],
                    dw_view,
                    threshold
                )
            end
        end
    end

    return nothing
end

# -------------------------------------------------------------------------
# 🚀 GENERIC WRAPPER
# -------------------------------------------------------------------------

"""
    calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, idx, type::Val, geom::Val; kwargs...)

Universal wrapper that dispatches the statistical calculation based on the
Temporal Type (Annual, Monthly, Seasonal) and Geometry Rank (Grid, Zonal, Global).

# Arguments
- `moa`, `soa`: Output buffers for Mean and Standard Deviation.
- `w_moa`: Output buffer for Aggregate Availability Weights.
- `ska`, `kua`: Placeholders for Skewness and Kurtosis.
- `data`: Input multi-dimensional tensor.
- `idx`: Indexing structure (Matrix or Array) corresponding to the `type`.
- `type::Val`: Temporal dispatch (e.g., `Val(:annual)`, `Val(:monthly)`, `Val(:seasonal)`).
- `geom::Val`: Geometry dispatch (e.g., `Val(:ij)`, `Val(:zonal)`, `Val(:global)`).
- `kwargs...`: Supports `weights`, `data_w`, and `threshold`.

# Returns
- `nothing` (Mutates output arrays in-place).
"""
calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, idx, type::Val, geom::Val; kwargs...) =
    (calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, idx, type; kwargs...); return nothing)

"""
    calc_anomalies!(anom_out, data, clim_mean, indices)

Calculates the departure from the climatological mean (anomalies) for each
time step in a simulation.

# Arguments
- `anom_out`: Output array for anomalies. Shape matches `data`.
- `data`: Original time series. Shape: (Sim, Var, [Spatial...], Time).
- `clim_mean`: Climatological mean (from Val(:monthly)).
   Shape: (Sim, Var, [Spatial...], 12).
- `indices`: Matrix{Vector{Int}} of size (Sim, 12).
   Mapping of which time-indices in `data` belong to which month (1-12).

# Logic
For each simulation and each month of the year (1-12), it identifies all
corresponding time-steps and performs: Anomaly = Data - Mean_of_that_month.
"""
function calc_anomalies!(anom_out, data, clim_mean, indices::AbstractMatrix)
    # Dimensionality check
    if ndims(anom_out) != ndims(data)
        throw(DimensionMismatch("Output and Input data must have same rank."))
    end

    Threads.@threads for ns in axes(data, 1)
        for m in 1:12
            t_vec = indices[ns, m]
            if !isempty(t_vec)
                # Slice the climatological mean for this specific month
                # (Sim, Var, [Spatial...], Month) -> (Var, [Spatial...])
                m_clim = selectdim(selectdim(clim_mean, 1, ns), ndims(clim_mean)-1, m)

                # Slice the data and output for all years of this month
                # (Sim, Var, [Spatial...], Time) -> (Var, [Spatial...], Selected_Years)
                d_slice = selectdim(selectdim(data, 1, ns), ndims(data)-1, t_vec)
                a_slice = selectdim(selectdim(anom_out, 1, ns), ndims(anom_out)-1, t_vec)

                # Broadcast the subtraction across the years
                a_slice .= d_slice .- m_clim
            end
        end
    end

    return nothing
end # function calc_anomalies!

"""
    calc_linear_trends!(slopes, data, weights, dates)

Calculates the weighted linear trend (slope) of a time series for every grid point.

# Arguments
- `slopes`: Output array for the trend. Shape: (Sim, Var, [Spatial...]).
- `data`: Input time series (typically anomalies). Shape: (Sim, Var, [Spatial...], Time).
- `weights`: Input weights (data_w). Shape: (Sim, Var, [Spatial...], Time).
- `dates`: Vector{Date} used to calculate the time axis (in years).

# Logic
Uses Weighted Least Squares (WLS). For each grid point, it solves for 'b' in:
    y = a + b*t
where weights 'w' account for data availability and month length.
"""
function calc_linear_trends!(slopes, data, weights, dates)
    t_axis = ndims(data)
    n_sims = size(data, 1)

    # Convert dates to fractional years (e.g., 2020.5) for a meaningful slope
    t_vals = Float32.([(year(d) + (month(d)-1)/12.0) for d in dates])

    Threads.@threads for ns in 1:n_sims
        # Select simulation slices
        d_sim = selectdim(data, 1, ns)
        w_sim = selectdim(weights, 1, ns)
        s_sim = selectdim(slopes, 1, ns)

        # Spatial dimensions loop (flattened for generic rank handling)
        # We iterate over all indices except the last (time) dimension
        for I in CartesianIndices(size(d_sim)[1:end-1])
            y = d_sim[I, :]
            w = w_sim[I, :]

            # Filter for valid data
            mask = .!isnan.(y) .& (w .> 0.0f0)
            if sum(mask) < 2
                s_sim[I] = NaN32
                continue
            end

            # Local Weighted Least Squares
            y_v = y[mask]
            t_v = t_vals[mask]
            w_v = w[mask]

            # WLS Formula for slope (b):
            # b = [ΣwΣwty - ΣwtΣwy] / [ΣwΣwt² - (Σwt)²]
            sw   = sum(w_v)
            swt  = sum(w_v .* t_v)
            swy  = sum(w_v .* y_v)
            swtt = sum(w_v .* t_v.^2)
            swty = sum(w_v .* t_v .* y_v)

            denom = (sw * swtt - swt^2)
            if denom ≈ 0
                s_sim[I] = NaN32
            else
                s_sim[I] = (sw * swty - swt * swy) / denom
            end
        end
    end

    return nothing
end # function calc_linear_trends!

"""
    calc_linear_trends_sig!(slopes, p_values, data, weights, dates)

Calculates the weighted linear trend and the statistical significance (p-value)
of that trend for every grid point.

# Arguments
- `slopes`: Output array for the slope (dT/dt).
- `p_values`: Output array for the p-values (0.0 to 1.0).
- `data`: Input anomalies (Sim, Var, [Spatial...], Time).
- `weights`: Input availability weights (Sim, Var, [Spatial...], Time).
- `dates`: Chronological Date vector.

# Logic
1. Solves Weighted Least Squares (WLS).
2. Calculates the Standard Error of the slope.
3. Performs a two-tailed t-test to determine if the slope is significantly different from zero.
"""
function calc_linear_trends_sig!(slopes, p_values, data, weights, dates)
    t_axis = ndims(data)
    n_sims = size(data, 1)

    # Time in fractional years
    t_vals = Float32.([(year(d) + (month(d)-1)/12.0) for d in dates])

    Threads.@threads for ns in 1:n_sims
        d_sim = selectdim(data, 1, ns)
        w_sim = selectdim(weights, 1, ns)
        s_sim = selectdim(slopes, 1, ns)
        p_sim = selectdim(p_values, 1, ns)

        for I in CartesianIndices(size(d_sim)[1:end-1])
            y = d_sim[I, :]
            w = w_sim[I, :]

            mask = .!isnan.(y) .& (w .> 0.0f0)
            n_valid = sum(mask)

            # We need at least 3 points to have a degree of freedom for error (n-2)
            if n_valid < 3
                s_sim[I], p_sim[I] = NaN32, NaN32
                continue
            end

            y_v, t_v, w_v = y[mask], t_vals[mask], w[mask]

            # WLS Components
            sw, swt, swy = sum(w_v), sum(w_v .* t_v), sum(w_v .* y_v)
            swtt, swty = sum(w_v .* t_v.^2), sum(w_v .* t_v .* y_v)

            delta = sw * swtt - swt^2
            if delta ≈ 0
                s_sim[I], p_sim[I] = NaN32, NaN32
                continue
            end

            # Slope (b) and Intercept (a)
            b = (sw * swty - swt * swy) / delta
            a = (swy - b * swt) / sw

            # Residual Variance Calculation
            y_pred = a .+ b .* t_v
            residuals = y_v .- y_pred
            # Weighted Residual Sum of Squares
            wrss = sum(w_v .* residuals.^2)

            # Standard Error of the slope
            # σ_b = sqrt( [WRSS / (n-2)] / Σ(w(t - t_avg)^2) )
            df = n_valid - 2
            mse = wrss / df
            # Variance of t: Σw(t - t_avg)² = Δ / Σw
            se_b = sqrt(mse / (delta / sw))

            # t-statistic
            t_stat = b / se_b

            # Store results
            s_sim[I] = b
            # Calculate p-value using a Student's T distribution
            # Note: For speed, one can use an approximation for large N,
            # but Distributions.jl is precise.
            tdist = TDist(df)
            p_sim[I] = 2 * ccdf(tdist, abs(t_stat))
        end
    end

    return nothing
end # function calc_linear_trends_sig!

"""
    generate_significance_mask!(mask, p_values; alpha=0.05)

Generates a binary mask where 1 indicates statistical significance and 0
indicates non-significance (noise).

# Arguments
- `mask`: Output array (Float32 or Int). Shape must match `p_values`.
- `p_values`: Input array of p-values (0.0 to 1.0) from trend analysis.
- `alpha`: Significance threshold (default 0.05).
  * p < alpha => Mask = 1 (Significant)
  * p >= alpha => Mask = 0 (Not Significant)

# Usage in Plotting
One can multiply the trend map by this mask to "hide" non-significant trends,
or use it to add stippling (dots) to the visualizations.
"""
function generate_significance_mask!(mask, p_values; alpha=0.05)
    # Check if p_values is empty or unallocated
    if isempty(p_values)
        @warn "p_values array is empty. No mask generated."
        return nothing
    end

    # Efficient broadcasting for the binary threshold
    # 1.0 represents 'True' (Significant), 0.0 represents 'False' (Noise)
    mask .= ifelse.(p_values .< Float32(alpha), 1.0f0, 0.0f0)

    # Handle NaNs: If the p-value is NaN (due to missing data), the mask is 0
    mask[isnan.(p_values)] .= 0.0f0

    return nothing
end # function generate_significance_mask!

"""
    calc_ensemble_stats!(ens_mean, ens_std, moa;
                         ens_ska=nothing, ens_kua=nothing,
                         w_moa=nothing, threshold=0.0)

Calculates weighted ensemble statistics across multiple simulations.

# Arguments
- `ens_mean`, `ens_std`: **Output** Ensemble Mean and Standard Deviation.
  Shape: `(Var, [Spatial...], [Time])`.
- `ens_ska`, `ens_kua`: **Output (Optional)** Ensemble Skewness and Kurtosis.
  Set to `nothing` if not needed.
- `moa`: **Input** Multi-year Means from individual simulations.
  Shape: `(Sim, Var, [Spatial...], [Time])`.
- `w_moa`: **Input (Optional)** Availability weights. Shape must match `moa`.
- `threshold`: Percentage (0-100) of cumulative availability required across
  the ensemble to produce a valid result.

# Mathematical Logic
1. If `w_moa` is provided, calculates weighted moments (Mean, Variance, Skew, Kurtosis).
2. If `w_moa` is `nothing`, performs standard unweighted arithmetic statistics.
3. Higher moments (Skewness \$G_1\$, Kurtosis \$G_2\$) are calculated using
   weighted standardized moments.
"""
function calc_ensemble_stats!(ens_mean, ens_std, moa;
                              ens_ska=nothing, ens_kua=nothing,
                              w_moa=nothing, threshold=0.0)
    n_sims = size(moa, 1)
    t_dim = ndims(moa)
    el_type = eltype(moa)

    # 1. Handle Optional Weights
    # If no weights provided, use 1.0 for non-NaN points
    weights = isnothing(w_moa) ? ifelse.(isnan.(moa), 0.0f0, 1.0f0) : w_moa
    sum_w = sum(weights, dims=1)

    # 2. Weighted Ensemble Mean (μ)
    mu_ens = _nansum(moa .* weights, dims=1) ./ sum_w

    # 3. Deviations and Variance (σ²)
    diff = moa .- mu_ens
    var_ens = _nansum(weights .* diff.^2, dims=1) ./ sum_w
    sigma_ens = sqrt.(var_ens)

    # 4. Optional Higher Moments
    # Skewness (3rd moment)
    if !isnothing(ens_ska)
        m3 = _nansum(weights .* diff.^3, dims=1) ./ sum_w
        ska = m3 ./ (sigma_ens.^3)
        # Threshold and collapse
        invalid_ska = (sum_w ./ Float32(n_sims) .* 100.0f0) .< threshold
        ens_ska .= dropdims(ifelse.(invalid_ska, el_type(NaN), ska), dims=1)
    end

    # Kurtosis (4th moment)
    if !isnothing(ens_kua)
        m4 = _nansum(weights .* diff.^4, dims=1) ./ sum_w
        kua = (m4 ./ (var_ens.^2)) .- 3.0f0 # Excess Kurtosis
        # Threshold and collapse
        invalid_kua = (sum_w ./ Float32(n_sims) .* 100.0f0) .< threshold
        ens_kua .= dropdims(ifelse.(invalid_kua, el_type(NaN), kua), dims=1)
    end

    # 5. Thresholding and Assignment for Mean/Std
    combined_fraction = sum_w ./ Float32(n_sims)
    invalid = combined_fraction .< (Float32(threshold) / 100.0f0)

    ens_mean .= dropdims(ifelse.(invalid, el_type(NaN), mu_ens), dims=1)
    ens_std .= dropdims(ifelse.(invalid, el_type(NaN), sigma_ens), dims=1)

    return nothing
end # function calc_ensemble_stats!

"""
    calc_ensemble_agreement!(pos_count, neg_count, agreement_frac, moa;
                             w_moa=nothing, threshold=0.0)

Calculates the level of consensus across an ensemble regarding the sign of
the data (e.g., positive vs. negative trends).

# Arguments
- `pos_count::AbstractArray`: **Output** Number of simulations with values > 0.
- `neg_count::AbstractArray`: **Output** Number of simulations with values < 0.
- `agreement_frac::AbstractArray`: **Output** The consensus strength [0.0 to 1.0].
  Calculated as: max(pos_count, neg_count) / total_valid_sims.
- `moa::AbstractArray`: **Input** Multi-year results (e.g., trends or anomalies).
  Shape: `(Sim, Var, [Spatial...], [Time])`.
- `w_moa::AbstractArray`: **Input (Optional)** Availability weights.
- `threshold::Float64`: Percentage (0-100) of total simulations required to
  have valid data for an agreement score to be calculated.

# Returns
- `nothing` (Mutates output arrays in-place).
"""
function calc_ensemble_agreement!(pos_count, neg_count, agreement_frac, moa;
                                  w_moa=nothing, threshold=0.0)
    n_sims = size(moa, 1)
    el_type = eltype(moa)

    # 1. Define valid mask (taking weights into account if provided)
    # A simulation is "valid" if it is not NaN and has a weight > 0
    valid_mask = isnothing(w_moa) ? .!isnan.(moa) : (w_moa .> 0.0f0) .& .!isnan.(moa)

    # 2. Count Positive and Negative signs across the 'Sim' dimension (dim 1)
    # We use sum(...) on booleans which Julia handles efficiently
    p_cnt = sum(valid_mask .& (moa .> 0), dims=1)
    n_cnt = sum(valid_mask .& (moa .< 0), dims=1)
    v_cnt = sum(valid_mask, dims=1) # Total valid simulations per point

    # 3. Calculate Agreement Fraction
    # If 9 out of 10 models are positive, agreement is 0.9.
    # If 5 are positive and 5 are negative, agreement is 0.5 (no consensus).
    max_agreement = max.(p_cnt, n_cnt)
    frac = max_agreement ./ v_cnt

    # 4. Thresholding
    # Mask out points where the number of contributing models is too low
    insufficient = (v_cnt ./ Float32(n_sims) .* 100.0f0) .< threshold

    # 5. Final Assignment
    pos_count .= dropdims(ifelse.(insufficient, 0, p_cnt), dims=1)
    neg_count .= dropdims(ifelse.(insufficient, 0, n_cnt), dims=1)
    agreement_frac .= dropdims(ifelse.(insufficient, el_type(NaN), frac), dims=1)

    return nothing
end # function calc_ensemble_agreement!

end # module PlanetStatistics
