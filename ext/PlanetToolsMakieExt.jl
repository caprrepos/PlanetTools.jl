## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-04 to 2026-10-31

# Package extension with the Makie-based methods of PlanetTools; loaded automatically
# once Makie (which brings ColorSchemes and Contour) is loaded together with PlanetTools.
module PlanetToolsMakieExt

using ColorSchemes: ColorScheme, colorschemes
using Contour: Contour
using Makie: Makie, cgrad, lines!
using PlanetTools.Constructors: Constructors
using PlanetTools.PlottingTools: PlottingTools

# --- Constructors.compose_palette ---

# Return a ColorScheme for a palette given as a ColorScheme, Symbol, or String.
_scheme(p::ColorScheme) = p
function _scheme(p::Union{Symbol,AbstractString})
    key = Symbol(p)
    if !haskey(colorschemes, key)
        msg = "unknown color scheme :$key (see keys(ColorSchemes.colorschemes))"
        throw(ArgumentError(msg))
    end
    return colorschemes[key]
end

# `n` colors sampled evenly over [0, 1] of a scheme; a single color is the scheme's
# midpoint (0.5) rather than its first color.
function _sample(p, n::Integer)
    n <= 0 && return Makie.RGBAf[]
    cs = _scheme(p)
    x = n == 1 ? (0.5,) : range(0.0, 1.0; length=n)
    return Makie.RGBAf[Makie.to_color(get(cs, xi)) for xi in x]
end

function Constructors.compose_palette(
    breaks::AbstractVector{<:Real};
    palette1=nothing,
    midcolor=nothing,
    palette2=nothing,
    demarcation::Union{Real,Nothing}=1.0,
)
    n = length(breaks)
    n >= 2 || throw(ArgumentError("compose_palette: need ≥ 2 breaks, got $n"))
    if !(all(isfinite, breaks) && all(i -> breaks[i] < breaks[i + 1], 1:(n - 1)))
        msg = "compose_palette: breaks must be finite and strictly increasing"
        throw(ArgumentError(msg))
    end
    if palette1 === nothing && palette2 === nothing
        msg = "compose_palette: at least one of palette1, palette2 must be given"
        throw(ArgumentError(msg))
    end

    d = demarcation === nothing ? 1.0 : demarcation
    nint = n - 1  # Number of colors to return.
    hasmid = midcolor !== nothing

    if palette1 !== nothing && palette2 !== nothing
        k = count(<=(d), breaks)
        nlow = max(k - 1, 0)
        strad = 1 <= k <= n - 1  # An interval crossing the demarcation exists.
        if hasmid && !strad
            @warn "compose_palette: no interval straddles demarcation = $d; " *
                "midcolor is not used"
        end
        nmid = hasmid && strad ? 1 : 0
        nhigh = nint - nlow - nmid
    elseif palette1 !== nothing
        nmid = hasmid ? 1 : 0
        nlow, nhigh = nint - nmid, 0
    else
        nmid = hasmid ? 1 : 0
        nlow, nhigh = 0, nint - nmid
    end

    pal = Makie.RGBAf[]
    sizehint!(pal, nint)
    append!(pal, _sample(palette1, nlow))
    nmid == 1 && push!(pal, Makie.to_color(midcolor))
    append!(pal, _sample(palette2, nhigh))

    @assert length(pal) == nint
    return pal
end

# --- PlottingTools.finite_contour! ---

function PlottingTools.finite_contour!(
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

function PlottingTools.finite_contour!(
    ax, x::AbstractVector, y::AbstractVector, z::AbstractMatrix{<:Real}, levels; kwargs...
)
    # The positional levels win over a levels keyword.
    return PlottingTools.finite_contour!(ax, x, y, z; kwargs..., levels=levels)
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

end # module PlanetToolsMakieExt
