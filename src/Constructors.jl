## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-11 to 2026-10-24

"""
    PlanetTools.Constructors

Submodule of PlanetTools with utility constructors for ModelE / ROCKE-3D post-processing
and plotting: [`fnames_modele`](@ref) builds monthly output file names and
[`compose_palette`](@ref) builds discrete color palettes for binned plots. Load it with
`using PlanetTools.Constructors`.
"""
module Constructors

using ColorSchemes: ColorScheme, colorschemes
using Dates: Date, Month, month, year
using Makie: Makie

export fnames_modele, compose_palette

# Locale-independent upper-case month abbreviations used in ModelE file names.
const MODELE_MONTHS = (
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
)

"""
    fnames_modele(fyear, fmonth, lyear, lmonth, gname, ids; ext=".nc") -> Matrix{String}

Build ModelE / ROCKE-3D monthly output file names `"<MON><YYYY>.<gname><id><ext>"`, e.g.
`"JAN1950.aijE001.nc"`, for each run ID in `ids` and each month from `fyear`-`fmonth` to
`lyear`-`lmonth` (inclusive). `<MON>` is the upper-case English month abbreviation
(locale-independent) and `<YYYY>` the year zero-padded to four digits.

# Arguments
- `fyear::Integer`: first year, `≥ 0`.
- `fmonth::Integer`: first month, `1:12`.
- `lyear::Integer`: last year, `≥ 0`.
- `lmonth::Integer`: last month, `1:12`.
- `gname::AbstractString`: diagnostic file group, e.g. `"aij"` or `"taij"`.
- `ids::AbstractVector{<:AbstractString}`: run IDs.

# Keywords
- `ext::AbstractString=".nc"`: file extension; `""` for none.

# Returns
- `Matrix{String}`: size `(length(ids), nmonths)`; rows are runs, columns are months in
    chronological order.

# Throws
- `ArgumentError`: if a month is outside `1:12`, a year is negative, the start date is
    after the end date, or `ids` is empty.

See also [`compose_palette`](@ref).

# Examples
```jldoctest
julia> fnames_modele(1950, 12, 1951, 1, "aij", ["E001", "E002"])
2×2 Matrix{String}:
 "DEC1950.aijE001.nc"  "JAN1951.aijE001.nc"
 "DEC1950.aijE002.nc"  "JAN1951.aijE002.nc"
```
"""
function fnames_modele(
    fyear::Integer,
    fmonth::Integer,
    lyear::Integer,
    lmonth::Integer,
    gname::AbstractString,
    ids::AbstractVector{<:AbstractString};
    ext::AbstractString=".nc",
)
    1 <= fmonth <= 12 || throw(ArgumentError("fmonth must be in 1:12, got $fmonth"))
    1 <= lmonth <= 12 || throw(ArgumentError("lmonth must be in 1:12, got $lmonth"))
    if fyear < 0 || lyear < 0
        throw(ArgumentError("years must be non-negative, got fyear=$fyear, lyear=$lyear"))
    end
    isempty(ids) && throw(ArgumentError("ids must contain at least one run ID"))

    dfirst = Date(fyear, fmonth)
    dlast = Date(lyear, lmonth)
    if dfirst > dlast
        msg = "start date $fyear-$fmonth is after end date $lyear-$lmonth"
        throw(ArgumentError(msg))
    end

    # "MONYYYY" stamps in chronological order.
    cmy = [MODELE_MONTHS[month(d)] * lpad(year(d), 4, '0') for d in dfirst:Month(1):dlast]

    # Rows are run IDs, columns are months.
    return [string(m, ".", gname, id, ext) for id in ids, m in cmy]
end

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

"""
    compose_palette(breaks; kwargs...) -> Vector{Makie.RGBAf}

Return one color per interval between consecutive `breaks` (`length(breaks) - 1` colors)
for discrete contour or binned plots, built from one or two color schemes and an optional
`midcolor`.

With both palettes, the intervals at or below `demarcation` take `palette1`, the interval
from the last break `≤ demarcation` to the next break takes `midcolor` (or else
`palette2`), and the intervals above take `palette2`. If that interval does not exist,
`midcolor` is not used and a warning is issued. With one palette, `midcolor` takes the
top interval (`palette1`) or the bottom interval (`palette2`). Each palette is sampled
evenly from its first to its last color; a single sample is its midpoint.

Use `contourf(x, y, z; levels=breaks, colormap=pal)`. For extension bands pass
`extendlow=first(pal), extendhigh=last(pal)`; `:auto` makes Makie resample the colormap.

# Arguments
- `breaks::AbstractVector{<:Real}`: finite, strictly increasing, at least two values.

# Keywords
- `palette1=nothing`: lower palette, a `ColorScheme` or its name as `Symbol`/`String`.
- `midcolor=nothing`: any color accepted by `Makie.to_color`, e.g. `:white`.
- `palette2=nothing`: upper palette, same types as `palette1`.
- `demarcation::Real=1.0`: value separating the palettes; used only if both are given.

# Returns
- `Vector{Makie.RGBAf}`: colors from the lowest to the highest interval.

# Throws
- `ArgumentError`: if `breaks` is invalid, no palette is given, or a scheme name is
    unknown.

See also [`fnames_modele`](@ref).

# Examples
```jldoctest
julia> pal = compose_palette(
           [-2, -1, -0.5, 0.5, 1, 2];
           palette1=reverse(colorschemes[:Blues]),
           midcolor=:white,
           palette2=:Reds,
           demarcation=0,
       );

julia> length(pal), pal[3] == Makie.to_color(:white)
(5, true)
```
"""
function compose_palette(
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

end # module Constructors
