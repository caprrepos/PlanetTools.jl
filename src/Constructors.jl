## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-11 to 2026-10-25

"""
    PlanetTools.Constructors

Submodule of PlanetTools with utility constructors for ModelE / ROCKE-3D post-processing
and plotting: [`fnames_modele`](@ref) builds monthly output file names and
[`compose_palette`](@ref) builds discrete color palettes for binned plots. Load it with
`using PlanetTools.Constructors`. `compose_palette` is defined in a package extension that
is loaded once Makie is loaded (`using Makie`, `using CairoMakie`, …).
"""
module Constructors

using Dates: Date, Month, month, year

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

"""
    compose_palette(breaks; kwargs...) -> Vector{Makie.RGBAf}

Return one color per interval between consecutive `breaks` (`length(breaks) - 1` colors)
for discrete contour or binned plots, built from one or two color schemes and an optional
`midcolor`. Requires Makie to be loaded.

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
function compose_palette end

# Hint for calls before the Makie extension is loaded.
function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        ext = Base.get_extension(parentmodule(@__MODULE__), :PlanetToolsMakieExt)
        if exc.f === compose_palette && ext === nothing
            print(
                io, "\n`compose_palette` is defined in the Makie extension of PlanetTools: "
            )
            print(io, "run `using Makie` (or CairoMakie, GLMakie) first.")
        end
    end
    return nothing
end

end # module Constructors
