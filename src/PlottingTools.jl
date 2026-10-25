## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-04 to 2026-10-25

"""
    PlanetTools.PlottingTools

Submodule of PlanetTools with plotting helpers for Makie (including GeoMakie axes):
[`finite_contour!`](@ref) draws isolines of fields with missing data. Load it with
`using PlanetTools.PlottingTools`. The methods are defined in a package extension that is
loaded once Makie is loaded (`using Makie`, `using CairoMakie`, …).
"""
module PlottingTools

export finite_contour!

"""
    finite_contour!(ax, x, y, z[, levels]; kwargs...) -> NamedTuple

Draw isolines of `z(x, y)` on the Makie axis `ax` (e.g. an `Axis` or a GeoMakie `GeoAxis`),
robustly for fields with missing data (`NaN`). Requires Makie to be loaded.

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
function finite_contour! end

# Hint for calls before the Makie extension is loaded.
function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        ext = Base.get_extension(parentmodule(@__MODULE__), :PlanetToolsMakieExt)
        if exc.f === finite_contour! && ext === nothing
            print(
                io, "\n`finite_contour!` is defined in the Makie extension of PlanetTools: "
            )
            print(io, "run `using Makie` (or CairoMakie, GLMakie) first.")
        end
    end
    return nothing
end

end # module PlottingTools
