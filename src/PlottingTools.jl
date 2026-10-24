## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-04 to 2026-10-24

"""
    PlanetTools.PlottingTools

Submodule of PlanetTools with plotting helpers for Makie (including GeoMakie axes):
[`finite_contour!`](@ref) draws isolines of fields with missing data. Load it with
`using PlanetTools.PlottingTools`.
"""
module PlottingTools

using Contour: Contour
using Makie: cgrad, lines!

export finite_contour!

"""
    finite_contour!(ax, x, y, z[, levels]; kwargs...) -> NamedTuple

Draw isolines of `z(x, y)` on the Makie axis `ax` (e.g. an `Axis` or a GeoMakie `GeoAxis`),
robustly for fields with missing data (`NaN`).

The isolines of each level are traced with Contour.jl (marching squares); vertices with
non-finite coordinates, which arise next to `NaN` cells, become gaps, and segments with
fewer than `minpoints` finite vertices are dropped. All segments of a level are drawn with
one `lines!` call, so the axis projection and line attributes apply as usual. Unlike
`Makie.contour!`, a segment without finite vertices can never cause an error.

# Arguments
- `ax`: target Makie axis.
- `x::AbstractVector`: coordinates of the first dimension of `z`, e.g. longitudes.
- `y::AbstractVector`: coordinates of the second dimension of `z`, e.g. latitudes.
- `z::AbstractMatrix{<:Real}`: field of size `(length(x), length(y))`; `NaN` marks missing
    cells.
- `levels`: as the keyword; a positional value takes precedence.

# Keywords
- `levels=5`: an `Integer` `n` gives `n` levels equally spaced strictly inside the finite
    range of `z` (or `zrange`); a vector gives these levels; any other `Real` gives one
    level. Note that `levels=1` is one automatic level, `levels=1.0` the isovalue 1.
- `color=nothing`: one color for all levels, a vector with one color per level (in the
    order of `levels`), or `nothing` to color the levels from `colormap`.
- `colormap=:viridis`: colormap for `color=nothing`, anything accepted by `cgrad`.
- `colorrange=nothing`: `(low, high)` mapped onto `colormap`; default the level extrema.
- `zrange=nothing`: `(low, high)` for automatic levels instead of the finite range of `z`,
    e.g. to exclude sentinel values.
- `minpoints::Integer=2`: minimum number of finite vertices of a drawn segment.
- `kwargs...`: passed to `lines!`, e.g. `linewidth` or `linestyle`.

# Returns
- `NamedTuple`: `levels`, the sorted unique levels, and `nsegments`, the number of segments
    drawn per level (`0` if a level does not occur). Automatic levels are empty if `z` has
    no finite values or a constant finite range.

# Throws
- `DimensionMismatch`: if `size(z) != (length(x), length(y))`.
- `ArgumentError`: for fewer than one automatic level, non-finite levels, an invalid
    `zrange`, `colorrange`, or `minpoints`, or a color vector of the wrong length.

# Examples
```jldoctest
julia> x, y = 0.0:4.0, 0.0:3.0;

julia> z = [xi for xi in x, yi in y];  # Isolines are vertical lines x = level.

julia> z[3, 3:4] .= NaN;  # Missing cells at x = 2, y ≥ 2.

julia> ax = Axis(Figure()[1, 1]);

julia> finite_contour!(ax, x, y, z, [2.5, 0.5]; color=[:red, :blue], linewidth=2)
(levels = [0.5, 2.5], nsegments = [1, 1])
```
"""
function finite_contour!(
    ax,
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractMatrix{<:Real};
    levels=5,
    color=nothing,
    colormap=:viridis,
    colorrange=nothing,
    zrange=nothing,
    minpoints::Integer=2,
    kwargs...,
)
    if size(z) != (length(x), length(y))
        nx, ny = length(x), length(y)
        msg = "size(z) = $(size(z)) but (length(x), length(y)) = ($nx, $ny)"
        throw(DimensionMismatch(msg))
    end
    minpoints >= 1 || throw(ArgumentError("minpoints must be ≥ 1, got $minpoints"))
    levs, order, nuser = _resolve_levels(levels, z, zrange)
    colors = _level_colors(color, colormap, colorrange, levs, order, nuser)
    nseg = zeros(Int, length(levs))

    for (il, lev) in pairs(levs)
        xs, ys = Float64[], Float64[]
        for line in Contour.lines(Contour.contour(x, y, z, lev))
            lx, ly = Contour.coordinates(line)
            good = isfinite.(lx) .& isfinite.(ly)
            count(good) < minpoints && continue
            append!(xs, ifelse.(good, lx, NaN), (NaN,))
            append!(ys, ifelse.(good, ly, NaN), (NaN,))
            nseg[il] += 1
        end
        nseg[il] > 0 && lines!(ax, xs, ys; color=colors[il], kwargs...)
    end
    return (; levels=levs, nsegments=nseg)
end

function finite_contour!(
    ax, x::AbstractVector, y::AbstractVector, z::AbstractMatrix{<:Real}, levels; kwargs...
)
    return finite_contour!(ax, x, y, z; kwargs..., levels=levels)  # Positional wins.
end

# Sorted unique levels, the index of each into the levels as given (for per-level colors),
# and the number of levels as given.
function _resolve_levels(n::Integer, z, zrange)
    n >= 1 || throw(ArgumentError("number of automatic levels must be ≥ 1, got $n"))
    lo, hi = zrange === nothing ? _finite_extrema(z) : _checked_range(zrange, "zrange")
    hi > lo || return Float64[], Int[], n  # Nothing to contour.
    return collect(range(lo, hi; length=n + 2))[2:(end - 1)], collect(1:n), n
end

function _resolve_levels(levels::AbstractVector{<:Real}, z, zrange)
    all(isfinite, levels) || throw(ArgumentError("levels must be finite"))
    order = sortperm(levels)
    unique!(i -> levels[i], order)  # Keep the first of equal levels.
    return Float64.(levels[order]), order, length(levels)
end

_resolve_levels(level::Real, z, zrange) = _resolve_levels([level], z, zrange)

function _finite_extrema(z)
    any(isfinite, z) || return (NaN, NaN)
    lo, hi = extrema(Iterators.filter(isfinite, z))
    return Float64(lo), Float64(hi)
end

function _checked_range(r, name)
    if !(length(r) == 2 && all(isfinite, r) && r[1] <= r[2])
        throw(ArgumentError("$name must be finite (low, high) with low ≤ high"))
    end
    return Float64(r[1]), Float64(r[2])
end

# One color per resolved level.
function _level_colors(color, colormap, colorrange, levs, order, nuser)
    if color === nothing
        isempty(levs) && return []
        cg = cgrad(colormap)
        lo, hi = if colorrange === nothing
            extrema(levs)
        else
            _checked_range(colorrange, "colorrange")
        end
        return [
            hi > lo ? get(cg, clamp((v - lo) / (hi - lo), 0, 1)) : get(cg, 0.5) for
            v in levs
        ]
    elseif color isa AbstractVector
        if length(color) != nuser
            throw(ArgumentError("got $(length(color)) colors for $nuser levels"))
        end
        return color[order]
    else
        return fill(color, length(levs))
    end
end

end # module PlottingTools
