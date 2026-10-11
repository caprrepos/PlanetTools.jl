"""
    Constructors

Small utility constructors used in ModelE / ROCKE-3D post-processing and plotting
workflows.

# Exported functions

| Function             | Purpose                                                                  |
|:---------------------|:-------------------------------------------------------------------------|
| [`fnames_modelE`](@ref) | Build the matrix of monthly ModelE output file names (`MONYYYY.<group><runID><ext>`) for several runs and a date range. |
| [`compose_palette`](@ref)  | Build a discrete color palette (one color per interval between contour/bin `breaks`) from one or two `ColorSchemes.jl` palettes and an optional middle color. |

# Dependencies

- `Dates` (standard library): month/year date ranges.
- `ColorSchemes.jl`: named color schemes (`colorschemes[:Blues]`, …).
- `Makie.jl`: only `Makie.to_color` and the `Makie.RGBAf` color type, so the
  returned palette can be passed directly to Makie plotting functions.

# Example

```julia
using .Constructors
fnames = fnames_modelE(2000, 1, 2000, 12, "aij", ["E6TomaF40", "E6TomaF41"])
pal    = compose_palette([-2, -1, -0.5, 0.5, 1, 2];
                         palette1 = :Blues, midcolor = :white,
                         palette2 = :Reds,  demarcation = 0.0)
```
# Notes

Anthropic's Claude Opus 5.5 Medium was used for debugging, code optimization, and writing docstrings
on Oct 11, 2026.
"""
module Constructors

export fnames_modelE, compose_palette

using ColorSchemes, Dates
import Makie

# Locale-independent upper-case month abbreviations used in ModelE file names.
const MODELE_MONTHS = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                       "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

"""
    fnames_modelE(fyear, fmonth, lyear, lmonth, gname, ids; ext = ".nc") -> Matrix{String}

Construct the names of monthly ModelE (or ROCKE-3D) post-processed output files for
every combination of simulation (run ID) and month in an inclusive date range.

Each file name is assembled as

    <MON><YYYY>.<gname><runID><ext>

for example `"JAN1950.aijE001.nc"`, where

- `<MON>`   is the upper-case three-letter English month abbreviation
  (`JAN`, `FEB`, …, `DEC`), independent of the system locale,
- `<YYYY>`  is the year zero-padded to (at least) four digits
  (year 850 → `"0850"`), matching the original implementation,
- `<gname>` is the diagnostic file group (e.g. `"aij"`, `"aijl"`, `"ijk"`,
  `"taij"`, `"adiurn"`),
- `<runID>` is the simulation identifier as it appears in the file name,
- `<ext>`   is the file extension (default `".nc"`).

The function only builds strings; it does not check whether the files exist on disk.
Prepend a directory with e.g. `joinpath.(datadir, fnames)`.

# Arguments
- `fyear::Integer`: first year of the range.
- `fmonth::Integer`: first month of the range, `1 ≤ fmonth ≤ 12`.
- `lyear::Integer`: last year of the range (inclusive).
- `lmonth::Integer`: last month of the range (inclusive), `1 ≤ lmonth ≤ 12`.
- `gname::AbstractString`: diagnostic file group name (e.g. `"aij"`).
- `ids::AbstractVector{<:AbstractString}`: run IDs, one row of the output per ID.
  Any string-like element type is accepted (e.g. the `SubString`s returned by
  `split`).

# Keywords
- `ext::AbstractString = ".nc"`: file extension appended to every name. Pass `""`
  to omit the extension.

# Returns
- `Matrix{String}` of size `(length(ids), nmonths)`: element `[i, j]` is the file
  for run `ids[i]` and the `j`-th month of the range, in chronological order.
  `nmonths = 12*(lyear - fyear) + (lmonth - fmonth) + 1`.

# Throws
- `ArgumentError` if a month is outside `1:12`.
- `ArgumentError` if the start date `(fyear, fmonth)` is later than the end date
  `(lyear, lmonth)` (the original version silently returned an empty matrix).
- `ArgumentError` if `ids` is empty.

# Example
```julia
julia> using .Constructors

julia> fnames_modelE(1950, 11, 1951, 2, "aij", ["E001", "E002"])
2×4 Matrix{String}:
 "NOV1950.aijE001.nc"  "DEC1950.aijE001.nc"  "JAN1951.aijE001.nc"  "FEB1951.aijE001.nc"
 "NOV1950.aijE002.nc"  "DEC1950.aijE002.nc"  "JAN1951.aijE002.nc"  "FEB1951.aijE002.nc"

julia> files = joinpath.("/discover/nobackup/me/E001", fnames_modelE(2000, 1, 2000, 12, "taij", ["E001"]; ext = ".nc"));

julia> vec(files)  # 12 monthly paths for one run
```
"""
function fnames_modelE(fyear::Integer, fmonth::Integer, lyear::Integer,
        lmonth::Integer, gname::AbstractString,
        ids::AbstractVector{<:AbstractString}; ext::AbstractString = ".nc")

    (1 <= fmonth <= 12) || throw(ArgumentError("fmonth must be in 1:12, got $fmonth"))
    (1 <= lmonth <= 12) || throw(ArgumentError("lmonth must be in 1:12, got $lmonth"))
    isempty(ids) && throw(ArgumentError("ids must contain at least one run ID"))

    dfirst = Date(fyear, fmonth)
    dlast  = Date(lyear, lmonth)
    dfirst <= dlast || throw(ArgumentError(
        "start date $(fyear)-$(fmonth) is after end date $(lyear)-$(lmonth)"))

    # "MONYYYY" stamps, chronological
    cmy = [MODELE_MONTHS[month(d)] * lpad(year(d), 4, '0')
           for d in dfirst:Month(1):dlast]

    # rows = run IDs, columns = months
    return [string(m, ".", gname, id, ext) for id in ids, m in cmy]
end # fnames_modelE


# Return a ColorScheme for a palette given as a ColorScheme, Symbol, or String.
_scheme(p::ColorScheme) = p
_scheme(p::Union{Symbol, AbstractString}) = colorschemes[Symbol(p)]

# n colors sampled evenly over [0, 1] of a scheme; a single color is the scheme's
# midpoint (0.5), as in the original intention of the "c == 1" branch.
function _sample(p, n::Integer)
    n <= 0 && return Makie.RGBAf[]
    cs = _scheme(p)
    x  = n == 1 ? [0.5] : collect(range(0.0, 1.0; length = n))
    return Makie.RGBAf[Makie.to_color(c) for c in get(cs, x)]
end

"""
    compose_palette(breaks; palette1 = nothing, midcolor = nothing,
                    palette2 = nothing, demarcation = nothing) -> Vector{Makie.RGBAf}

Compose a discrete color palette with exactly one color per interval between
consecutive `breaks` (i.e. `length(breaks) - 1` colors), from one or two continuous
color schemes and an optional middle color. Typical use: discrete filled contours or
binned heat maps (`contourf`, `heatmap` with categorical bins), especially diverging
plots of anomalies (`demarcation = 0`) or ratios (`demarcation = 1`, the default).

# How the colors are assigned

Let `n = length(breaks)`, so there are `n - 1` intervals
`[breaks[i], breaks[i+1]]`, `i = 1, …, n-1`.

**Two palettes (`palette1` and `palette2`)** — diverging palette. With
`k = count(≤(demarcation), breaks)`:

- the `max(k - 1, 0)` intervals lying entirely at or below `demarcation` are colored
  with `palette1` (sampled from its start, value 0, to its end, value 1);
- the interval that straddles `demarcation` (from the last break `≤ demarcation`
  to the first break `> demarcation`; it exists only if `1 ≤ k ≤ n-1`) gets
  `midcolor` if one is given, otherwise it is the first `palette2` interval;
- the remaining intervals above get `palette2` (sampled from 0 to 1).

To center `midcolor` on the demarcation value, choose `breaks` so that the
demarcation value lies *strictly between* two breaks, e.g.
`breaks = [-2, -1, -0.5, 0.5, 1, 2]` with `demarcation = 0`. If `demarcation`
equals a break, the "straddling" interval is the one just above it.
If no break lies on one side of `demarcation`, there is no straddling interval:
all intervals go to the other palette and `midcolor` is not used (a warning is
issued).

**Only `palette1`** — `midcolor` (if given) becomes the *last* color (top interval);
the other `n - 2` (or `n - 1` without `midcolor`) colors come from `palette1`.

**Only `palette2`** — `midcolor` (if given) becomes the *first* color (bottom
interval); the rest come from `palette2`.

A palette that contributes exactly one color contributes its midpoint color
(`get(scheme, 0.5)`). Each palette is always sampled from its first to its last
color; pass a reversed scheme (`reverse(colorschemes[:Blues])`) to flip a branch.
For a conventional blue–white–red anomaly map, `palette1` is usually a *reversed*
blue scheme so that the darkest blue sits in the lowest interval.

# Arguments
- `breaks::AbstractVector{<:Real}`: interval boundaries, sorted in non-decreasing
  order, `length(breaks) ≥ 2`.

# Keywords
- `palette1 = nothing`: lower palette; a `Symbol`/`String` naming a scheme in
  `ColorSchemes.colorschemes` (e.g. `:Blues`, `:viridis`) or a `ColorScheme`.
- `midcolor = nothing`: any color accepted by `Makie.to_color`
  (`:white`, `"gray90"`, `RGBf(1,1,1)`, …).
- `palette2 = nothing`: upper palette, same accepted types as `palette1`.
- `demarcation = nothing`: value separating the two palettes; `nothing` means `1.0`
  (useful for ratios). Only used when both palettes are given.

At least one of `palette1`, `palette2` must be given.

# Returns
- `Vector{Makie.RGBAf}` of length `length(breaks) - 1`, ordered from the lowest to
  the highest interval — directly usable as `colormap` (with
  `Makie.Categorical` / `cgrad(pal; categorical = true)`) or for `contourf` with
  `levels = breaks`.

# Throws
- `ArgumentError` if `length(breaks) < 2`, `breaks` is not sorted, or no palette is
  given.

# Example
```julia
using .Constructors, CairoMakie, ColorSchemes

# Diverging anomaly palette: 2 blues | white around 0 | 2 reds  (5 intervals)
breaks = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]
pal = compose_palette(breaks;
                      palette1 = reverse(colorschemes[:Blues]),
                      midcolor = :white,
                      palette2 = :Reds,
                      demarcation = 0.0)
length(pal) == length(breaks) - 1        # true

x = range(-3, 3; length = 200)
z = [sin(a) * cos(b) * 2 for a in x, b in x]
fig, ax, cf = contourf(x, x, z; levels = breaks, colormap = pal,
                       extendlow = :auto, extendhigh = :auto)
Colorbar(fig[1, 2], cf)
fig

# Sequential palette with a white top class
compose_palette(0:10:100; palette1 = :viridis, midcolor = :white)   # 10 colors
```
"""
function compose_palette(breaks; palette1 = nothing, midcolor = nothing,
        palette2 = nothing, demarcation = nothing)

    n = length(breaks)
    n >= 2 || throw(ArgumentError(
        "compose_palette: number of breaks must be ≥ 2, got $n"))
    issorted(breaks) || throw(ArgumentError(
        "compose_palette: breaks must be sorted in non-decreasing order"))
    (palette1 === nothing && palette2 === nothing) && throw(ArgumentError(
        "compose_palette: at least one of palette1, palette2 must be given"))

    d    = demarcation === nothing ? 1.0 : demarcation
    nint = n - 1                       # number of colors to return
    hasmid = midcolor !== nothing

    if palette1 !== nothing && palette2 !== nothing
        k     = count(<=(d), breaks)
        nlow  = max(k - 1, 0)
        strad = 1 <= k <= n - 1        # interval crossing the demarcation exists
        if hasmid && !strad
            @warn "compose_palette: no interval straddles demarcation = $d; midcolor is not used"
        end
        nmid  = (hasmid && strad) ? 1 : 0
        nhigh = nint - nlow - nmid
    elseif palette1 !== nothing
        nmid  = hasmid ? 1 : 0
        nlow, nhigh = nint - nmid, 0
    else
        nmid  = hasmid ? 1 : 0
        nlow, nhigh = 0, nint - nmid
    end

    pal = Makie.RGBAf[]
    append!(pal, _sample(palette1, nlow))
    nmid == 1 && push!(pal, Makie.to_color(midcolor))
    append!(pal, _sample(palette2, nhigh))

    @assert length(pal) == nint
    return pal
end # compose_palette

end # module Constructors
