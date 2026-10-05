module PlottingTools

export finite_contour!

using Contour                       # Contour.jl (marching squares)
using Makie: lines!, cgrad          # also works with GLMakie/CairoMakie + GeoMakie

# ───────────────────────── internal helpers ─────────────────────────
### created with Anthropic's Claude Opus 5.5, October 4, 2026

"""
    _finite_extrema(z) -> (zmin::Float64, zmax::Float64)

Return the minimum and maximum over the *finite* entries of `z`, ignoring
`NaN` and `±Inf`. Throws `ArgumentError` if `z` has no finite entry.

# Example
```julia
julia> _finite_extrema([NaN 1.0; Inf -2.0])
(-2.0, 1.0)
```
"""
function _finite_extrema(z)
    any(isfinite, z) ||
        throw(ArgumentError("finite_contour!: z contains no finite values"))
    lo, hi = extrema(Iterators.filter(isfinite, z))
    return Float64(lo), Float64(hi)
end

"""
    _resolve_levels(levels, z, zrange) -> Vector{Float64}

Turn the user's `levels` specification into an explicit, sorted, duplicate-free
vector of isovalues.

- `levels::Integer` (n ≥ 1): `n` equally spaced levels strictly inside
  `(zmin, zmax)`, i.e. `range(zmin, zmax; length = n + 2)[2:end-1]`. Here
  `(zmin, zmax)` is `zrange` if given, otherwise the finite extrema of `z`.
- `levels::AbstractVector{<:Real}`: these exact levels.
- `levels::AbstractFloat`: this one level.

# Example
```julia
julia> _resolve_levels(3, [0.0 4.0; 1.0 2.0], nothing)
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```
"""
function _resolve_levels(n::Integer, z, zrange)
    n ≥ 1 || throw(ArgumentError("finite_contour!: number of levels must be ≥ 1, got $n"))
    lo, hi = zrange === nothing ? _finite_extrema(z) : (Float64(zrange[1]), Float64(zrange[2]))
    hi > lo || throw(ArgumentError("finite_contour!: degenerate range ($lo, $hi) for automatic levels"))
    return collect(range(lo, hi; length = n + 2))[2:end-1]
end
_resolve_levels(levels::AbstractVector{<:Real}, z, zrange) = sort!(unique(Float64.(levels)))
_resolve_levels(level::AbstractFloat, z, zrange) = [Float64(level)]

"""
    _level_colors(color, colormap, colorrange, levs) -> Vector

Return one colour per level in `levs`.

- `color === nothing`: sample `colormap` at each level's position within
  `colorrange` (default: `extrema(levs)`; a single level gets the middle
  of the colormap).
- `color::AbstractVector`: used as is; its length must equal `length(levs)`.
- anything else (`:yellow`, `"red"`, `RGBf(…)`, `(:black, 0.5)`): repeated
  for every level.

# Example
```julia
julia> length(_level_colors(nothing, :viridis, nothing, [1.0, 2.0, 3.0]))
3
```
"""
function _level_colors(color, colormap, colorrange, levs)
    nl = length(levs)
    if color === nothing
        cg = cgrad(colormap)
        lo, hi = colorrange === nothing ? extrema(levs) : Float64.(colorrange)
        return [hi > lo ? get(cg, clamp((v - lo) / (hi - lo), 0, 1)) : get(cg, 0.5) for v in levs]
    elseif color isa AbstractVector
        length(color) == nl || throw(ArgumentError(
            "finite_contour!: got $(length(color)) colors for $nl levels"))
        return collect(color)
    else
        return fill(color, nl)
    end
end

# ───────────────────────────── public API ────────────────────────────
### created with Anthropic's Claude Opus 5.5, October 4, 2026

"""
    finite_contour!(ax, x, y, z; levels = 5, color = nothing, colormap = :viridis,
                    colorrange = nothing, zrange = nothing, minpoints = 2, kwargs...)
        -> (; levels::Vector{Float64}, nsegments::Vector{Int})

    finite_contour!(ax, x, y, z, levels; kwargs...)      # positional form

Draw isolines of `z(x, y)` on the Makie axis `ax` (for example a GeoMakie
`GeoAxis`) in a way that stays robust when `z` contains missing data
(`NaN`). It works like `Makie.contour!`, but it cannot fail with
`ArgumentError: reducing over an empty collection is not allowed`, an
error Makie's `contour!` can raise when non-finite values produce isoline
segments without finite vertices.

# Algorithm
For every requested level:
1. trace the isolines with `Contour.contour(x, y, z, level)` (marching
   squares);
2. mark vertices with non-finite coordinates as `NaN` (pen-up gaps);
3. drop segments with fewer than `minpoints` finite vertices;
4. join the remaining segments, separated by `NaN`, into one coordinate
   vector, and draw it with a single `lines!` call, so that the axis' map
   projection and all line attributes apply as usual.

# Positional arguments
- `ax`: target Makie axis (`Axis`, `GeoAxis`, …).
- `x::AbstractVector`: coordinates of the first dimension of `z`
  (longitudes of the cell centres), `length(x) == size(z, 1)`.
- `y::AbstractVector`: coordinates of the second dimension of `z`
  (latitudes), `length(y) == size(z, 2)`.
- `z::AbstractMatrix{<:Real}`: field to contour. `NaN` marks undefined
  cells. Replace `±Inf` by large finite values beforehand if those cells
  carry an "above or below the level" meaning.
- `levels` (optional fifth positional argument): same meaning as the
  keyword. When given, it takes precedence over the keyword.

# Keywords
- `levels = 5`: the level specification, chosen by type:
  * `Integer n` → `n` automatic levels, equally spaced strictly inside
    the finite range of `z` (or `zrange`);
  * `AbstractVector{<:Real}` → exactly these levels (sorted, duplicates
    removed);
  * `AbstractFloat` → this one level.
  Note: `levels = 1` means *one automatic level*; write `levels = 1.0` or
  `levels = [1]` for the isovalue 1.
- `color = nothing`: a single colour for all levels, a vector with one
  colour per level, or `nothing` to colour each level from `colormap`.
- `colormap = :viridis`: any colormap accepted by `cgrad`; used only if
  `color === nothing`.
- `colorrange = nothing`: `(lo, hi)` value range mapped onto `colormap`;
  defaults to the extrema of the resolved levels.
- `zrange = nothing`: `(zmin, zmax)` used for automatic levels instead of
  the finite extrema of `z`. Useful when sentinel values (e.g.
  `floatmax(Float32)` replacing `Inf`) would otherwise distort the range.
- `minpoints::Int = 2`: minimum number of finite vertices for a segment to
  be drawn.
- `kwargs...`: passed on to `lines!` (`linewidth`, `linestyle`, `alpha`,
  `label`, …).

# Returns
A `NamedTuple` `(; levels, nsegments)` giving the resolved levels and the
number of segments drawn for each. `nsegments[i] == 0` means level `i`
does not occur in the finite part of `z`.

# Throws
- `DimensionMismatch` if `size(z) != (length(x), length(y))`;
- `ArgumentError` for a non-positive number of levels, a degenerate
  automatic range, an all-non-finite `z`, or a colour vector whose length
  does not match the number of levels.

# Example
```julia
using GLMakie, GeoMakie, Contour

lon = Float32.(-177.5:5:177.5)
lat = Float32.(-90:4:90)
z   = Float32[cosd(φ) * (1 + sind(2λ)) for λ in lon, φ in lat]
z[10:15, 20:25] .= NaN32                         # undefined cells

fig = Figure(size = (900, 900))

# (a) one explicit level, fixed colour
ga1 = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = "single level")
finite_contour!(ga1, lon, lat, z; levels = [0.5], color = :yellow, linewidth = 0.5)

# (b) several explicit levels, one colour each
ga2 = GeoAxis(fig[2, 1]; dest = "+proj=moll", title = "explicit levels")
finite_contour!(ga2, lon, lat, z; levels = [0.25, 0.5, 1.0],
                color = [:blue, :black, :red], linewidth = 1)

# (c) 8 automatic levels coloured by value
ga3 = GeoAxis(fig[3, 1]; dest = "+proj=moll", title = "8 automatic levels")
res = finite_contour!(ga3, lon, lat, z; levels = 8, colormap = :plasma)
@info "levels drawn" res.levels res.nsegments
fig
```
Bibliography

Source type: all entries are primary package documentation. None of them has a DOI.

[1] JuliaGeometry. Contour.jl documentation (contour, contours, levels, lines, coordinates). No DOI. https://juliageometry.github.io/Contour.jl/stable/

[2] MakieOrg. Makie documentation: contour (an Int value for levels means a number of levels, a vector means explicit levels). No DOI. https://docs.makie.org/stable/reference/plots/contour

[3] MakieOrg. Makie documentation: lines (NaN separates segments) and Colors / cgrad. No DOI. https://docs.makie.org/stable/reference/plots/lines ; https://docs.makie.org/stable/explanations/colors
"""
function finite_contour!(ax, x::AbstractVector, y::AbstractVector, z::AbstractMatrix{<:Real};
                         levels = 5, color = nothing, colormap = :viridis,
                         colorrange = nothing, zrange = nothing,
                         minpoints::Int = 2, kwargs...)
    size(z) == (length(x), length(y)) ||
        throw(DimensionMismatch("size(z) = $(size(z)) but (length(x), length(y)) = " *
                                "($(length(x)), $(length(y)))"))
    levs   = _resolve_levels(levels, z, zrange)
    colors = _level_colors(color, colormap, colorrange, levs)
    nseg   = zeros(Int, length(levs))

    for (il, lev) in pairs(levs)
        xs_all = Float32[]; ys_all = Float32[]
        for ln in Contour.lines(Contour.contour(x, y, z, lev))
            xs, ys = Contour.coordinates(ln)
            good = isfinite.(xs) .& isfinite.(ys)
            count(good) < minpoints && continue
            append!(xs_all, ifelse.(good, Float32.(xs), NaN32), (NaN32,))
            append!(ys_all, ifelse.(good, Float32.(ys), NaN32), (NaN32,))
            nseg[il] += 1
        end
        nseg[il] > 0 && lines!(ax, xs_all, ys_all; color = colors[il], kwargs...)
    end
    return (; levels = levs, nsegments = nseg)
end

finite_contour!(ax, x::AbstractVector, y::AbstractVector, z::AbstractMatrix{<:Real},
                levels; kwargs...) =
    finite_contour!(ax, x, y, z; levels = levels, kwargs...)
end # module Plotting Tools
