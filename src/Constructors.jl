## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-11 to 2026-10-31

"""
    PlanetTools.Constructors

Submodule of PlanetTools with name constructors: [`fnames_modele`](@ref) builds monthly
output file names for ModelE / ROCKE-3D post-processing. Load it with
`using PlanetTools.Constructors`. `compose_palette`, formerly in this submodule, is now in
[`PlottingTools`](@ref PlanetTools.PlottingTools); its old name here is deprecated.
"""
module Constructors

using Dates: Date, Month, month, year
using ..PlottingTools: PlottingTools

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

# Moved to PlottingTools after v0.2.1; the old name keeps working and is deprecated.
# Remove with the next breaking release.
Base.@deprecate_binding compose_palette PlottingTools.compose_palette

end # module Constructors
