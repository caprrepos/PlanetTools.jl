## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-26

"""
    PlanetTools.PlanetStatistics

Submodule of PlanetTools with statistics for gridded planetary model output: zonal and
global means, annual and seasonal mean time series, multi-year statistics, anomalies,
linear trends with significance, and ensemble statistics. Load it with
`using PlanetTools.PlanetStatistics`.

Conventions shared by all functions:
- Arrays are laid out as `(Sim, Var, Lon, Lat, [Lev], Time)`; reduced dimensions are kept
    as singletons (e.g. zonal means `(Sim, Var, 1, Lat, [Lev], Time)`).
- Missing data are `NaN`.
- Optional availability weights (fractions in `[0, 1]`, e.g. returned by an earlier step)
    multiply the temporal weights (e.g. days per month); points with weight `0` count as
    missing.
- `threshold` is the minimum percentage of valid points for a result; otherwise the result
    is `NaN` and its availability weight `0`.
- Results are written into preallocated output arrays; entries not covered by a
    calculation are left unchanged, so initialize outputs with `NaN`.
- Simulations are processed in parallel with `Threads.@threads`; within a simulation the
    output points are processed in parallel as well, so a single simulation also uses all
    threads.
- The kernels are explicit loops that accumulate in `Float64` and allocate no temporary
    arrays of the size of the data.
- Standard deviations are sample standard deviations (`corrected=true`, Bessel's
    correction generalized to weights) unless `corrected=false` is passed.
"""
module PlanetStatistics

using Dates: Date, month, year
using SpecialFunctions: beta_inc

export calc_planet_means!, calc_annual_mean_ts!, calc_seasonal_mean_ts!
export calc_multiyear_stats!, calc_anomalies!
export calc_linear_trends!, calc_linear_trends_sig!, generate_significance_mask!
export calc_ensemble_stats!, calc_ensemble_agreement!

const MONTHS_PER_YEAR = 12

"""
    calc_planet_means!(
        data, zm_out, gm_out, weights_zm, weights_gm, axyp, t_ranges; kwargs...
    ) -> Nothing

Compute zonal means and area-weighted global means of `data` and the availability of valid
data for each.

The zonal mean averages over longitude (dimension 3) with equal weights; its availability
is the fraction of valid longitudes. The global mean averages over longitude and latitude
(dimensions 3 and 4) weighted by the cell areas `axyp`; its availability is the fraction of
the area with valid data.

# Arguments
- `data::AbstractArray`: `(Sim, Var, Lon, Lat, [Lev], Time)`.
- `zm_out`, `weights_zm`: outputs, `(Sim, Var, 1, Lat, [Lev], Time)`.
- `gm_out`, `weights_gm`: outputs, `(Sim, Var, 1, 1, [Lev], Time)`.
- `axyp::AbstractMatrix`: cell areas, `(Lon, Lat)`.
- `t_ranges::AbstractVector{<:Integer}`: number of time steps of each simulation; steps
    `1:t_ranges[ns]` are processed.

# Keywords
- `threshold::Real=100.0`: minimum valid percentage (of longitudes or area).
- `run_names=nothing`: simulation names for log messages.
- `verbose::Bool=false`: log progress with `@info`.

# Examples
```jldoctest
julia> data = reshape([1.0, 3.0, 2.0, NaN], 1, 1, 2, 2, 1);  # 2 lon × 2 lat

julia> axyp = [1.0 3.0; 1.0 3.0];  # Second latitude has 3× larger cells.

julia> zm, wz = fill(NaN, 1, 1, 1, 2, 1), fill(NaN, 1, 1, 1, 2, 1);

julia> gm, wg = fill(NaN, 1, 1, 1, 1, 1), fill(NaN, 1, 1, 1, 1, 1);

julia> calc_planet_means!(data, zm, gm, wz, wg, axyp, [1]; threshold=0);

julia> vec(zm), vec(wz), gm[1], wg[1]
([2.0, 2.0], [1.0, 0.5], 2.0, 0.625)
```
"""
function calc_planet_means!(
    data::AbstractArray,
    zm_out::AbstractArray,
    gm_out::AbstractArray,
    weights_zm::AbstractArray,
    weights_gm::AbstractArray,
    axyp::AbstractMatrix,
    t_ranges::AbstractVector{<:Integer};
    threshold::Real=100.0,
    run_names=nothing,
    verbose::Bool=false,
)
    nd = ndims(data)
    nd >= 5 || throw(ArgumentError("data must have at least 5 dimensions"))
    if size(axyp) != (size(data, 3), size(data, 4))
        throw(DimensionMismatch("axyp must have size (Lon, Lat) of data"))
    end
    total_area = sum(axyp)
    Threads.@threads for ns in eachindex(t_ranges)
        _log(verbose, "Zonal and global means", ns, run_names)
        idx = (ns, ntuple(_ -> Colon(), nd - 2)..., 1:t_ranges[ns])
        _planet_means_kernel!(
            view(zm_out, idx...),
            view(weights_zm, idx...),
            view(gm_out, idx...),
            view(weights_gm, idx...),
            view(data, idx...),
            axyp,
            total_area,
            threshold,
        )
    end
    return nothing
end

"""
    calc_annual_mean_ts!(
        data, zm, gm, ann_out, zm_ann_out, gm_ann_out, ann_w_out, zm_ann_w_out,
        gm_ann_w_out, dates, t_ranges, full_year_dates; kwargs...
    ) -> Nothing

Compute weighted annual means of monthly grid-point, zonal-mean, and global-mean time
series, with the fraction of valid months per year as availability.

Year `ny` of simulation `ns` averages the 12 months starting at the first date of
`full_year_dates[ns]` plus `ny - 1` years; years extending beyond `t_ranges[ns]` are
skipped.
Monthly weights are the availability weights times `month_weights`.

# Arguments
- `data`, `zm`, `gm`: monthly inputs, `(Sim, Var, Lon, Lat, [Lev], Time)` and the zonal and
    global reductions thereof.
- `ann_out`, `zm_ann_out`, `gm_ann_out`: annual means, input shapes with `Time` replaced by
    `Year`.
- `ann_w_out`, `zm_ann_w_out`, `gm_ann_w_out`: annual availability, same shapes.
- `dates::AbstractVector{<:Date}`: dates of the `Time` dimension.
- `t_ranges::AbstractVector{<:Integer}`: number of valid time steps of each simulation.
- `full_year_dates::AbstractVector`: per simulation, a date range of the complete years to
    average (e.g. `Date(1951):Year(1):Date(2000)`); empty to skip a simulation.

# Keywords
- `data_w=nothing`, `zm_w=nothing`, `gm_w=nothing`: availability weights of the inputs.
- `month_weights=nothing`: weights along `Time`, e.g. days per month; default equal.
- `threshold::Real=100.0`: minimum percentage of valid months per year.
- `run_names=nothing`, `verbose::Bool=false`: progress logging.
"""
function calc_annual_mean_ts!(
    data::AbstractArray,
    zm::AbstractArray,
    gm::AbstractArray,
    ann_out::AbstractArray,
    zm_ann_out::AbstractArray,
    gm_ann_out::AbstractArray,
    ann_w_out::AbstractArray,
    zm_ann_w_out::AbstractArray,
    gm_ann_w_out::AbstractArray,
    dates::AbstractVector{<:Date},
    t_ranges::AbstractVector{<:Integer},
    full_year_dates::AbstractVector;
    data_w=nothing,
    zm_w=nothing,
    gm_w=nothing,
    month_weights=nothing,
    threshold::Real=100.0,
    run_names=nothing,
    verbose::Bool=false,
)
    fields = (
        (data, data_w, ann_out, ann_w_out),
        (zm, zm_w, zm_ann_out, zm_ann_w_out),
        (gm, gm_w, gm_ann_out, gm_ann_w_out),
    )
    for (x, xw, _, _) in fields
        _check_same_size(x, xw)
    end

    Threads.@threads for ns in eachindex(t_ranges)
        _log(verbose, "Annual means", ns, run_names)
        years = full_year_dates[ns]
        isempty(years) && continue
        start = first(years)
        i0 = findfirst(d -> year(d) == year(start) && month(d) == month(start), dates)
        if i0 === nothing
            @warn "Start date $start of simulation $ns not found in dates"
            continue
        end
        for ny in 1:length(years)
            t = (i0 + MONTHS_PER_YEAR * (ny - 1)):(i0 + MONTHS_PER_YEAR * ny - 1)
            last(t) > t_ranges[ns] && break
            tw = month_weights === nothing ? nothing : view(month_weights, t)
            for (x, xw, out, wout) in fields
                _time_mean!(
                    _slice(out, ns, ny),
                    _slice(wout, ns, ny),
                    _slice(x, ns, t),
                    xw === nothing ? nothing : _slice(xw, ns, t),
                    tw,
                    MONTHS_PER_YEAR,
                    threshold,
                )
            end
        end
    end
    return nothing
end

"""
    calc_seasonal_mean_ts!(
        data, zm, gm, out_ij, out_zm, out_gm, w_out_ij, w_out_zm, w_out_gm,
        ind_seas_first, dates; kwargs...
    ) -> Nothing

Compute weighted seasonal means of monthly grid-point, zonal-mean, and global-mean time
series, with the fraction of valid months per season as availability.

Season `ne` of year `ny` of simulation `ns` averages the `length_of_season` time steps
starting at `ind_seas_first[ns, ne, ny]`; an index of `0` ends that season's years.

# Arguments
- `data`, `zm`, `gm`: monthly inputs, `(Sim, Var, Lon, Lat, [Lev], Time)` and the zonal and
    global reductions thereof.
- `out_ij`, `out_zm`, `out_gm`: seasonal means, `(Sim, Var, Season, …, Year)` with the
    spatial dimensions of the inputs.
- `w_out_ij`, `w_out_zm`, `w_out_gm`: seasonal availability, same shapes.
- `ind_seas_first::AbstractArray{<:Integer,3}`: `(Sim, Season, Year)` index of the first
    month of each season in `Time`.
- `dates::AbstractVector{<:Date}`: dates of the `Time` dimension.

# Keywords
- `data_w=nothing`, `zm_w=nothing`, `gm_w=nothing`: availability weights of the inputs.
- `month_weights=nothing`: weights along `Time`, e.g. days per month; default equal.
- `length_of_season::Integer=3`: months per season.
- `threshold::Real=100.0`: minimum percentage of valid months per season.
- `run_names=nothing`, `verbose::Bool=false`: progress logging.
"""
function calc_seasonal_mean_ts!(
    data::AbstractArray,
    zm::AbstractArray,
    gm::AbstractArray,
    out_ij::AbstractArray,
    out_zm::AbstractArray,
    out_gm::AbstractArray,
    w_out_ij::AbstractArray,
    w_out_zm::AbstractArray,
    w_out_gm::AbstractArray,
    ind_seas_first::AbstractArray{<:Integer,3},
    dates::AbstractVector{<:Date};
    data_w=nothing,
    zm_w=nothing,
    gm_w=nothing,
    month_weights=nothing,
    length_of_season::Integer=3,
    threshold::Real=100.0,
    run_names=nothing,
    verbose::Bool=false,
)
    fields = (
        (data, data_w, out_ij, w_out_ij),
        (zm, zm_w, out_zm, w_out_zm),
        (gm, gm_w, out_gm, w_out_gm),
    )
    for (x, xw, _, _) in fields
        _check_same_size(x, xw)
    end
    if size(data, ndims(data)) != length(dates)
        throw(DimensionMismatch("length(dates) must equal the Time dimension of data"))
    end

    Threads.@threads for ns in axes(data, 1)
        _log(verbose, "Seasonal means", ns, run_names)
        for ne in axes(ind_seas_first, 2)
            for ny in axes(ind_seas_first, 3)
                t0 = ind_seas_first[ns, ne, ny]
                t0 == 0 && break  # No further years for this season.
                t = t0:(t0 + length_of_season - 1)
                tw = month_weights === nothing ? nothing : view(month_weights, t)
                for (x, xw, out, wout) in fields
                    _time_mean!(
                        _slice(out, ns, ne, ny),
                        _slice(wout, ns, ne, ny),
                        _slice(x, ns, t),
                        xw === nothing ? nothing : _slice(xw, ns, t),
                        tw,
                        length_of_season,
                        threshold,
                    )
                end
            end
        end
    end
    return nothing
end

"""
    calc_multiyear_stats!(moa, soa, w_moa, ska, kua, data, indices, kind::Val; kwargs...)

Compute weighted multi-year mean, standard deviation, and optionally skewness and excess
kurtosis over the time steps selected by `indices`, with the fraction of valid steps as
availability.

With weights `w` (availability × temporal weight) of the valid steps, `V₁ = Σw`,
`V₂ = Σw²`, the mean is `μ = Σwx / V₁` and the variance

- `corrected=true` (default): `σ² = Σw(x - μ)² / (V₁ - V₂ / V₁)`, the unbiased estimator
    for reliability weights; with equal weights this is `Σ(x - μ)² / (n - 1)`, the sample
    variance. It needs at least two valid steps (otherwise `NaN`).
- `corrected=false`: the population variance `Σw(x - μ)² / V₁`.

Skewness `m₃ / m₂^{3/2}` and excess kurtosis `m₄ / m₂² - 3` are always computed from the
population central moments `mₖ = Σw(x - μ)ᵏ / V₁`.

- `Val(:annual)`: `data` is `(Sim, Var, …, Time)`, outputs are `(Sim, Var, …)`, and
    `indices` is a `(Sim, 2)` matrix of the first and last time index.
- `Val(:monthly)`: `data` is `(Sim, Var, …, Time)`, outputs are `(Sim, Var, …, 12)`, and
    `indices` is a `(Sim, 12)` matrix of time-index vectors, one per calendar month.
- `Val(:seasonal)`: `data` is `(Sim, Var, Season, …, Year)`, outputs are
    `(Sim, Var, Season, …)`, and `indices` is a `(Sim, Season, 2)` array of the first and
    last year.

Simulations (or months, seasons) with empty or zero indices are skipped.

# Arguments
- `moa`, `soa`, `w_moa`: outputs for mean, standard deviation, and availability.
- `ska`, `kua`: outputs for skewness and excess kurtosis, or `nothing`.
- `data::AbstractArray`: input, see above.
- `indices`: selected time steps, see above.
- `kind::Val`: `Val(:annual)`, `Val(:monthly)`, or `Val(:seasonal)`.

# Keywords
- `weights=nothing`: `(Sim, Time)` (or `(Sim, Year)`) matrix of temporal weights.
- `data_w=nothing`: availability weights, same size as `data`.
- `threshold::Real=100.0`: minimum percentage of valid time steps.
- `corrected::Bool=true`: sample (`true`) or population (`false`) standard deviation.

# Examples
```jldoctest
julia> data = reshape([1.0, 2.0, 3.0, NaN], 1, 1, 4);  # (Sim, Var, Time)

julia> moa, soa, w = fill(NaN, 1, 1), fill(NaN, 1, 1), fill(NaN, 1, 1);

julia> calc_multiyear_stats!(moa, soa, w, nothing, nothing, data, [1 4], Val(:annual);
           threshold=50);

julia> moa[1], soa[1], w[1]  # Sample standard deviation of [1, 2, 3].
(2.0, 1.0, 0.75)

julia> calc_multiyear_stats!(moa, soa, w, nothing, nothing, data, [1 4], Val(:annual);
           threshold=50, corrected=false);

julia> soa[1] ≈ sqrt(2 / 3)  # Population standard deviation.
true
```
"""
function calc_multiyear_stats!(
    moa::AbstractArray,
    soa::AbstractArray,
    w_moa::AbstractArray,
    ska,
    kua,
    data::AbstractArray,
    indices::AbstractMatrix{<:Integer},
    ::Val{:annual};
    weights=nothing,
    data_w=nothing,
    threshold::Real=100.0,
    corrected::Bool=true,
)
    _check_same_size(data, data_w)
    Threads.@threads for ns in axes(data, 1)
        i1, i2 = indices[ns, 1], indices[ns, 2]
        (i1 > 0 && i2 >= i1) || continue
        _multiyear!(
            _slice.((moa, soa, w_moa, ska, kua), ns)...,
            _slice(data, ns),
            _slice(data_w, ns),
            i1:i2,
            weights === nothing ? nothing : view(weights, ns, i1:i2),
            threshold,
            corrected,
        )
    end
    return nothing
end

function calc_multiyear_stats!(
    moa::AbstractArray,
    soa::AbstractArray,
    w_moa::AbstractArray,
    ska,
    kua,
    data::AbstractArray,
    indices::AbstractMatrix{<:AbstractVector},
    ::Val{:monthly};
    weights=nothing,
    data_w=nothing,
    threshold::Real=100.0,
    corrected::Bool=true,
)
    _check_same_size(data, data_w)
    Threads.@threads for ns in axes(data, 1)
        for m in 1:MONTHS_PER_YEAR
            t = indices[ns, m]
            isempty(t) && continue
            outs = map(
                a -> a === nothing ? nothing : _lastdim(_slice(a, ns), m),
                (moa, soa, w_moa, ska, kua),
            )
            _multiyear!(
                outs...,
                _slice(data, ns),
                _slice(data_w, ns),
                t,
                weights === nothing ? nothing : view(weights, ns, t),
                threshold,
                corrected,
            )
        end
    end
    return nothing
end

function calc_multiyear_stats!(
    moa::AbstractArray,
    soa::AbstractArray,
    w_moa::AbstractArray,
    ska,
    kua,
    data::AbstractArray,
    indices::AbstractArray{<:Integer,3},
    ::Val{:seasonal};
    weights=nothing,
    data_w=nothing,
    threshold::Real=100.0,
    corrected::Bool=true,
)
    _check_same_size(data, data_w)
    Threads.@threads for ns in axes(data, 1)
        for s in axes(indices, 2)
            i1, i2 = indices[ns, s, 1], indices[ns, s, 2]
            (i1 > 0 && i2 >= i1) || continue
            season(a) = a === nothing ? nothing : selectdim(_slice(a, ns), 2, s)
            _multiyear!(
                season.((moa, soa, w_moa, ska, kua))...,
                season(data),
                season(data_w),
                i1:i2,
                weights === nothing ? nothing : view(weights, ns, i1:i2),
                threshold,
                corrected,
            )
        end
    end
    return nothing
end

"""
    calc_anomalies!(anom_out, data, clim_mean, indices) -> Nothing

Compute anomalies from a monthly climatology: for each simulation `ns` and month `m`,
`anom_out[ns, …, t] = data[ns, …, t] - clim_mean[ns, …, m]` for all `t in indices[ns, m]`.
Time steps not listed in `indices` are left unchanged.

# Arguments
- `anom_out`: output, same size as `data`.
- `data::AbstractArray`: `(Sim, Var, …, Time)`.
- `clim_mean::AbstractArray`: `(Sim, Var, …, 12)`, e.g. from `calc_multiyear_stats!` with
    `Val(:monthly)`.
- `indices::AbstractMatrix{<:AbstractVector}`: `(Sim, 12)` time indices of each month.
"""
function calc_anomalies!(
    anom_out::AbstractArray,
    data::AbstractArray,
    clim_mean::AbstractArray,
    indices::AbstractMatrix{<:AbstractVector},
)
    size(anom_out) == size(data) ||
        throw(DimensionMismatch("anom_out and data must have the same size"))
    if size(clim_mean)[1:(end - 1)] != size(data)[1:(end - 1)]
        throw(DimensionMismatch("clim_mean must be (Sim, Var, …, 12) matching data"))
    end
    Threads.@threads for ns in axes(data, 1)
        for m in 1:MONTHS_PER_YEAR
            t = indices[ns, m]
            isempty(t) && continue
            clim = _lastdim(_slice(clim_mean, ns), m)
            _lastdim(_slice(anom_out, ns), t) .= _lastdim(_slice(data, ns), t) .- clim
        end
    end
    return nothing
end

"""
    calc_linear_trends!(slopes, data, weights, dates) -> Nothing

Compute the weighted least-squares linear trend of each time series in `data`, in units of
`data` per year. Time is the fractional year `year + (month - 1) / 12` of `dates`.

At least two valid points (non-`NaN`, weight `> 0`) are required; otherwise the slope is
`NaN`. The fit is done in `Float64` with centered time.

# Arguments
- `slopes`: output, `(Sim, Var, …)`.
- `data::AbstractArray`: `(Sim, Var, …, Time)`, typically anomalies.
- `weights`: weights, same size as `data`, or `nothing` for equal weights.
- `dates::AbstractVector{<:Date}`: dates of the `Time` dimension.

See also [`calc_linear_trends_sig!`](@ref).

# Examples
```jldoctest
julia> dates = Date(2000, 1):Month(1):Date(2000, 12);

julia> data = reshape([0.5 * (i - 1) / 12 for i in 1:12], 1, 1, 12);

julia> slopes = fill(NaN, 1, 1);

julia> calc_linear_trends!(slopes, data, nothing, dates);

julia> slopes[1] ≈ 0.5
true
```
"""
function calc_linear_trends!(slopes::AbstractArray, data::AbstractArray, weights, dates)
    tyears = _check_trend_inputs(data, weights, dates)
    Threads.@threads for ns in axes(data, 1)
        x, w, s = _slice(data, ns), _slice(weights, ns), _slice(slopes, ns)
        for I in CartesianIndices(size(x)[1:(end - 1)])
            t, y, wv = _valid_series(tyears, x, w, I)
            s[I] = length(y) < 2 ? NaN : _wls(t, y, wv).slope
        end
    end
    return nothing
end

"""
    calc_linear_trends_sig!(slopes, p_values, data, weights, dates; kwargs...) -> Nothing

Compute weighted least-squares linear trends as [`calc_linear_trends!`](@ref) and their
two-sided p-values from a Student's t-test of zero slope with `n - 2` degrees of freedom.

With `autocorrelation=true`, the sample size is reduced to the effective sample size
`n_e = n (1 - r₁) / (1 + r₁)` (with `n_e ≤ n`), where `r₁` is the lag-1 autocorrelation of
the residuals in time order, following Santer et al. [2000]. Without it, autocorrelated
series such as monthly anomalies get p-values that are too small.

At least three valid points are required (`n_e > 2` with `autocorrelation=true`); otherwise
slope and p-value are `NaN`.

# Arguments
- `slopes`, `p_values`: outputs, `(Sim, Var, …)`.
- `data::AbstractArray`: `(Sim, Var, …, Time)`.
- `weights`: weights, same size as `data`, or `nothing` for equal weights.
- `dates::AbstractVector{<:Date}`: dates of the `Time` dimension.

# Keywords
- `autocorrelation::Bool=false`: correct for lag-1 autocorrelation of the residuals.

# References
- Santer, B. D., et al. (2000), Statistical significance of trends and trend differences in
    layer-average atmospheric temperature time series, *J. Geophys. Res.*, 105(D6),
    7337–7356, doi:10.1029/1999JD901105.

# Examples
```jldoctest
julia> dates = Date(2000, 1):Month(1):Date(2001, 8);

julia> data = reshape([0.1 * i + (-1)^i for i in 1:20], 1, 1, 20);

julia> slopes, p = fill(NaN, 1, 1), fill(NaN, 1, 1);

julia> calc_linear_trends_sig!(slopes, p, data, nothing, dates);

julia> round(slopes[1]; digits=3), round(p[1]; digits=4)
(1.38, 0.0112)
```
"""
function calc_linear_trends_sig!(
    slopes::AbstractArray,
    p_values::AbstractArray,
    data::AbstractArray,
    weights,
    dates;
    autocorrelation::Bool=false,
)
    tyears = _check_trend_inputs(data, weights, dates)
    Threads.@threads for ns in axes(data, 1)
        x, w = _slice(data, ns), _slice(weights, ns)
        s, p = _slice(slopes, ns), _slice(p_values, ns)
        for I in CartesianIndices(size(x)[1:(end - 1)])
            t, y, wv = _valid_series(tyears, x, w, I)
            s[I], p[I] = _trend_and_p(t, y, wv, autocorrelation)
        end
    end
    return nothing
end

"""
    generate_significance_mask!(mask, p_values; alpha=0.05) -> Nothing

Set `mask` to `1` where `p_values < alpha` (significant) and to `0` elsewhere, including
where `p_values` is `NaN`. The mask can be multiplied with a trend map or used for
stippling.

# Examples
```jldoctest
julia> mask = zeros(3);

julia> generate_significance_mask!(mask, [0.01, 0.2, NaN]);

julia> mask
3-element Vector{Float64}:
 1.0
 0.0
 0.0
```
"""
function generate_significance_mask!(
    mask::AbstractArray, p_values::AbstractArray; alpha::Real=0.05
)
    size(mask) == size(p_values) ||
        throw(DimensionMismatch("mask and p_values must have the same size"))
    mask .= ifelse.(p_values .< alpha, 1, 0)
    return nothing
end

"""
    calc_ensemble_stats!(ens_mean, ens_std, moa; kwargs...) -> Nothing

Compute the weighted ensemble mean and standard deviation over the simulations (dimension 1
of `moa`), and optionally skewness and excess kurtosis. The standard deviation is the
sample standard deviation for reliability weights (`corrected=true`, see
[`calc_multiyear_stats!`](@ref)) or the population standard deviation
(`corrected=false`); skewness and kurtosis are population moments.

The ensemble availability is `Σ w / n_sims` over valid members, with `w` from `w_moa`
(or 1); where it is below `threshold` percent, the results are `NaN`.

# Arguments
- `ens_mean`, `ens_std`: outputs, `moa` without dimension 1.
- `moa::AbstractArray`: `(Sim, Var, …)`, e.g. multi-year means of each simulation.

# Keywords
- `ens_ska=nothing`, `ens_kua=nothing`: outputs for skewness and excess kurtosis.
- `w_moa=nothing`: availability weights, same size as `moa`.
- `threshold::Real=0.0`: minimum ensemble availability in percent.
- `corrected::Bool=true`: sample (`true`) or population (`false`) standard deviation.

# Examples
```jldoctest
julia> moa = reshape([1.0, 2.0, 6.0], 3, 1);  # 3 simulations, 1 variable

julia> m, s = fill(NaN, 1), fill(NaN, 1);

julia> calc_ensemble_stats!(m, s, moa);

julia> m[1], s[1] ≈ sqrt(14 / 2)
(3.0, true)
```
"""
function calc_ensemble_stats!(
    ens_mean::AbstractArray,
    ens_std::AbstractArray,
    moa::AbstractArray;
    ens_ska=nothing,
    ens_kua=nothing,
    w_moa=nothing,
    threshold::Real=0.0,
    corrected::Bool=true,
)
    _check_same_size(moa, w_moa)
    higher = ens_ska !== nothing || ens_kua !== nothing
    st = _weighted_moments(moa, w_moa, nothing, 1, size(moa, 1); higher=higher)
    frac = st.sumw ./ size(moa, 1)
    invalid = dropdims(frac .* 100 .< threshold; dims=1)
    nanif(x) = ifelse.(invalid, NaN, dropdims(x; dims=1))
    ens_mean .= nanif(st.mean)
    var = corrected ? _corrected_var.(st.m2, st.sw, st.sw2) : st.m2
    ens_std .= nanif(sqrt.(var))
    ens_ska === nothing || (ens_ska .= nanif(st.m3 ./ st.m2 .^ 1.5))
    ens_kua === nothing || (ens_kua .= nanif(st.m4 ./ st.m2 .^ 2 .- 3))
    return nothing
end

"""
    calc_ensemble_agreement!(
        pos_count, neg_count, agreement_frac, moa; kwargs...
    ) -> Nothing

Count the simulations (dimension 1 of `moa`) with positive and negative values and compute
the sign agreement `max(pos, neg) / n_valid`, from 0.5 (no consensus) to 1 (all agree).
Valid members are non-`NaN` with weight `> 0`. Where fewer than `threshold` percent of the
simulations are valid, counts are 0 and the agreement is `NaN`.

# Arguments
- `pos_count`, `neg_count`, `agreement_frac`: outputs, `moa` without dimension 1.
- `moa::AbstractArray`: `(Sim, Var, …)`, e.g. trends of each simulation.

# Keywords
- `w_moa=nothing`: availability weights, same size as `moa`.
- `threshold::Real=0.0`: minimum percentage of valid simulations.

# Examples
```jldoctest
julia> moa = reshape([0.3, 0.1, -0.2, 0.4], 4, 1);

julia> p, n, a = zeros(Int, 1), zeros(Int, 1), fill(NaN, 1);

julia> calc_ensemble_agreement!(p, n, a, moa);

julia> p[1], n[1], a[1]
(3, 1, 0.75)
```
"""
function calc_ensemble_agreement!(
    pos_count::AbstractArray,
    neg_count::AbstractArray,
    agreement_frac::AbstractArray,
    moa::AbstractArray;
    w_moa=nothing,
    threshold::Real=0.0,
)
    _check_same_size(moa, w_moa)
    valid = w_moa === nothing ? .!isnan.(moa) : .!isnan.(moa) .& (w_moa .> 0)
    npos = sum(valid .& (moa .> 0); dims=1)
    nneg = sum(valid .& (moa .< 0); dims=1)
    nvalid = sum(valid; dims=1)
    insufficient = nvalid ./ size(moa, 1) .* 100 .< threshold
    pos_count .= dropdims(ifelse.(insufficient, 0, npos); dims=1)
    neg_count .= dropdims(ifelse.(insufficient, 0, nneg); dims=1)
    frac = ifelse.(insufficient, NaN, max.(npos, nneg) ./ nvalid)
    agreement_frac .= dropdims(frac; dims=1)
    return nothing
end

# --- Internal helpers ---

# Weighted moments of `x` along `dim` (used by the ensemble statistics). Valid points are
# non-NaN with positive availability `avail` (array like `x` or `nothing`); their weight is
# availability × temporal weight `tw` (vector along `dim` or `nothing`). Returns
# singleton-`dim` arrays: mean, population central moments m2 (and m3, m4 if `higher`), the
# valid fraction `count / nexpected`, the sum of availability weights `sumw` (temporal
# weights excluded), and the sums of the weights `sw` and of their squares `sw2`.
function _weighted_moments(x, avail, tw, dim, nexpected; higher::Bool=false)
    a = avail === nothing ? true : avail
    valid = .!isnan.(x) .& (a .> 0)
    twb = tw === nothing ? true : reshape(tw, ntuple(i -> i == dim ? length(tw) : 1, dim))
    w = ifelse.(valid, a .* twb, false)
    xv = ifelse.(valid, x, zero(eltype(x)))
    sw = sum(w; dims=dim)
    sw2 = sum(w .^ 2; dims=dim)
    mu = sum(w .* xv; dims=dim) ./ sw
    fraction = sum(valid; dims=dim) ./ nexpected
    sumw = sum(ifelse.(valid, a, false); dims=dim)
    d = ifelse.(valid, x .- mu, zero(eltype(mu)))
    m2 = sum(w .* d .^ 2; dims=dim) ./ sw
    m3 = m4 = nothing
    if higher
        m3 = sum(w .* d .^ 3; dims=dim) ./ sw
        m4 = sum(w .* d .^ 4; dims=dim) ./ sw
    end
    return (; mean=mu, m2, m3, m4, fraction, sumw, sw, sw2)
end

# Variance for reliability weights from the population variance `m2`, `V₁ = sw`, and
# `V₂ = sw2`: Σw(x - μ)² / (V₁ - V₂ / V₁); NaN if the denominator is not positive (fewer
# than two valid points).
function _corrected_var(m2, sw, sw2)
    denom = sw - sw2 / sw
    return denom > 0 ? m2 * sw / denom : oftype(float(m2), NaN)
end

# Apply `f(I)` to every index `I` of the `CartesianIndices` `R`, in parallel chunks when
# `R` is large enough to make threading worthwhile.
function _tforeach(f, R::CartesianIndices)
    n = length(R)
    nchunks = min(cld(n, 256), 4 * Threads.nthreads())
    if nchunks <= 1
        for I in R
            f(I)
        end
        return nothing
    end
    chunks = collect(Iterators.partition(1:n, cld(n, nchunks)))
    Threads.@threads for c in chunks
        for k in c
            f(R[k])
        end
    end
    return nothing
end

# Availability of `x[I, t]`: 1 without availability weights.
@inline _avail(::Nothing, I, t) = 1.0
@inline _avail(xw, I, t) = Float64(xw[I, t])

# Temporal weight of the k-th selected step: 1 without temporal weights.
@inline _tweight(::Nothing, k) = 1.0
@inline _tweight(tw, k) = Float64(tw[k])

# Store a value and its availability at `I`, applying the threshold (percent).
@inline function _store1!(out, wout, I, value, fraction, threshold)
    if fraction * 100 < threshold
        out[I] = NaN
        wout[I] = 0
    else
        out[I] = value
        wout[I] = fraction
    end
    return nothing
end

# Zonal (equal weights over longitude) and area-weighted global means of
# `x = (Var, Lon, Lat, rest...)` into `z`, `wz` (`(Var, 1, Lat, rest...)`) and `g`, `wg`
# (`(Var, 1, 1, rest...)`).
function _planet_means_kernel!(z, wz, g, wg, x, axyp, total_area, threshold)
    nv, nlon, nlat = size(x, 1), size(x, 2), size(x, 3)
    nhoriz = nlon * nlat
    _tforeach(CartesianIndices(size(x)[4:end])) do r
        @inbounds for v in 1:nv
            gs = 0.0
            ga = 0.0
            gn = 0
            for j in 1:nlat
                zs = 0.0
                zn = 0
                for i in 1:nlon
                    xi = x[v, i, j, r]
                    isnan(xi) && continue
                    a = Float64(axyp[i, j])
                    zs += xi
                    zn += 1
                    gs += a * xi
                    ga += a
                end
                zval = zn > 0 ? zs / zn : NaN
                _store1!(z, wz, CartesianIndex(v, 1, j, r), zval, zn / nlon, threshold)
                gn += zn
            end
            gfrac = gn == nhoriz ? 1.0 : ga / total_area  # Exactly 1 if complete.
            _store1!(g, wg, CartesianIndex(v, 1, 1, r), gs / ga, gfrac, threshold)
        end
    end
    return nothing
end

# Weighted mean of `x` over its last (time) dimension into `out`, availability into `wout`.
function _time_mean!(out, wout, x, xw, tw, nexpected, threshold)
    nt = size(x, ndims(x))
    _tforeach(CartesianIndices(out)) do I
        s = 0.0
        sw = 0.0
        n = 0
        @inbounds for k in 1:nt
            v = x[I, k]
            a = _avail(xw, I, k)
            (isnan(v) || !(a > 0)) && continue
            w = a * _tweight(tw, k)
            s += w * v
            sw += w
            n += 1
        end
        return _store1!(out, wout, I, s / sw, n / nexpected, threshold)
    end
    return nothing
end

# Multi-year statistics of `x[…, t]` over the last dimension of `x`.
function _multiyear!(moa, soa, w_moa, ska, kua, x, xw, t, tw, threshold, corrected)
    nt = length(t)
    higher = ska !== nothing || kua !== nothing
    _tforeach(CartesianIndices(moa)) do I
        s = 0.0
        sw = 0.0
        sw2 = 0.0
        n = 0
        @inbounds for (k, tk) in enumerate(t)
            v = x[I, tk]
            a = _avail(xw, I, tk)
            (isnan(v) || !(a > 0)) && continue
            w = a * _tweight(tw, k)
            s += w * v
            sw += w
            sw2 += w^2
            n += 1
        end
        fraction = n / nt
        if fraction * 100 < threshold
            moa[I] = NaN
            soa[I] = NaN
            w_moa[I] = 0
            ska === nothing || (ska[I] = NaN)
            kua === nothing || (kua[I] = NaN)
            return nothing
        end
        mu = s / sw
        m2 = 0.0
        m3 = 0.0
        m4 = 0.0
        @inbounds for (k, tk) in enumerate(t)
            v = x[I, tk]
            a = _avail(xw, I, tk)
            (isnan(v) || !(a > 0)) && continue
            w = a * _tweight(tw, k)
            d = v - mu
            m2 += w * d^2
            if higher
                m3 += w * d^3
                m4 += w * d^4
            end
        end
        var_pop = m2 / sw
        moa[I] = mu
        soa[I] = sqrt(corrected ? _corrected_var(var_pop, sw, sw2) : var_pop)
        w_moa[I] = fraction
        ska === nothing || (ska[I] = (m3 / sw) / var_pop^1.5)
        kua === nothing || (kua[I] = (m4 / sw) / var_pop^2 - 3)
        return nothing
    end
    return nothing
end

# Views of simulation `ns`: `_slice(a, ns)` drops dimension 1; `_slice(a, ns, j)` also
# selects `j` in the last dimension; `_slice(a, ns, s, j)` selects season `s` (dimension 3)
# and `j` in the last dimension.
_slice(::Nothing, args...) = nothing
_slice(a::AbstractArray, ns) = selectdim(a, 1, ns)
_slice(a::AbstractArray, ns, j) = view(a, ns, ntuple(_ -> Colon(), ndims(a) - 2)..., j)
function _slice(a::AbstractArray, ns, s, j)
    return view(a, ns, :, s, ntuple(_ -> Colon(), ndims(a) - 4)..., j)
end

# Select index (or indices) `j` of the last dimension.
_lastdim(a, j) = selectdim(a, ndims(a), j)

function _check_same_size(x, w)
    if w !== nothing && size(w) != size(x)
        throw(DimensionMismatch("weights of size $(size(w)) do not match data $(size(x))"))
    end
    return nothing
end

function _log(verbose, what, ns, run_names)
    if verbose
        name = run_names === nothing ? "simulation $ns" : run_names[ns]
        @info "$what: $name" thread = Threads.threadid()
    end
    return nothing
end

# --- Trends ---

function _check_trend_inputs(data, weights, dates)
    _check_same_size(data, weights)
    if length(dates) != size(data, ndims(data))
        throw(DimensionMismatch("length(dates) must equal the Time dimension of data"))
    end
    return [year(d) + (month(d) - 1) / MONTHS_PER_YEAR for d in dates]
end

# Valid (non-NaN, weight > 0) points of the series at `I`, as Float64 vectors.
function _valid_series(tyears, x, w, I)
    y = Float64.(x[I, :])
    wv = w === nothing ? ones(length(y)) : Float64.(w[I, :])
    keep = .!isnan.(y) .& (wv .> 0)
    return tyears[keep], y[keep], wv[keep]
end

# Weighted least squares with centered time.
function _wls(t, y, w)
    sw = sum(w)
    tm = sum(w .* t) / sw
    ym = sum(w .* y) / sw
    stt = sum(w .* (t .- tm) .^ 2)
    slope = stt > 0 ? sum(w .* (t .- tm) .* (y .- ym)) / stt : NaN
    return (; slope, intercept=ym - slope * tm, stt)
end

function _trend_and_p(t, y, w, autocorrelation)
    n = length(y)
    n < 3 && return NaN, NaN
    fit = _wls(t, y, w)
    isnan(fit.slope) && return NaN, NaN
    r = y .- (fit.intercept .+ fit.slope .* t)
    neff = float(n)
    if autocorrelation
        r1 = sum(r[1:(end - 1)] .* r[2:end]) / sum(r .^ 2)
        neff = isfinite(r1) ? min(n * (1 - r1) / (1 + r1), n) : neff
        neff > 2 || return fit.slope, NaN
    end
    df = neff - 2
    se = sqrt(sum(w .* r .^ 2) / df / fit.stt)
    p = se > 0 ? _two_sided_t_pvalue(fit.slope / se, df) : (fit.slope == 0 ? NaN : 0.0)
    return fit.slope, p
end

# Two-sided p-value of Student's t: P(|T| ≥ |t|) = I_{ν/(ν+t²)}(ν/2, 1/2).
_two_sided_t_pvalue(tstat, df) = first(beta_inc(df / 2, 0.5, df / (df + tstat^2)))

end # module PlanetStatistics
