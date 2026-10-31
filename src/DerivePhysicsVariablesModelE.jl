## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

"""
    PlanetTools.DerivePhysicsVariablesModelE

Submodule of PlanetTools that derives new physics variables from existing ones in ModelE /
ROCKE-3D output: sums, differences, and quotients of variables, and cumulative dust size
classes. Load it with `using PlanetTools.DerivePhysicsVariablesModelE`.

All functions use the same data layout: `sdata` is an array whose dimension 1 indexes the
simulations and dimension 2 the variables named in `vnames`; any further dimensions
(e.g. longitude, latitude, level, time) are carried along unchanged. `attrib`, if given,
holds one dictionary of NetCDF attributes per variable. The functions do not modify their
inputs; they return new `vnames`, `sdata`, and `attrib` with the derived variables
appended along dimension 2. Fill values should be converted to `missing` or `NaN` before.

Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
(jan.p.perlwitz@caprllc.com).
"""
module DerivePhysicsVariablesModelE

export accumulate_variables, subtract_variables, divide_variables
export accumulate_dust_masses_modele, accumulate_dust_optics_modele

# Dimension of `sdata` that indexes the variables.
const VDIM = 2

# Dust size bins and the cumulative classes built from them, in accumulation order.
const MASS_BINS = ("Clay", "Silt1", "Silt2", "Silt3", "Silt4", "Silt5")
const MASS_CLASSES = ("PM2Dust", "PM4Dust", "PM8Dust", "PM16Dust", "PM32Dust", "PM64Dust")
const OPTICS_BINS = (
    "Clay1", "Clay2", "Clay3", "Clay4", "Silt1", "Silt2", "Silt3", "Silt4", "Silt5"
)
const OPTICS_CLASSES = (
    "SumClay1Dust",
    "SumClay2Dust",
    "SumClay3Dust",
    "PM2Dust",
    "PM4Dust",
    "PM8Dust",
    "PM16Dust",
    "PM32Dust",
    "PM64Dust",
)

"""
    accumulate_variables(vnames, accnames, addname, sdata; kwargs...)

Return `vnames`, `sdata`, and `attrib` with a new variable `addname` appended: the sum of
all variables selected by `accnames`.

A `String` in `accnames` selects the variable with exactly that name; a `Regex` selects all
variables whose names match it. A variable selected more than once is summed once.

# Arguments
- `vnames::AbstractVector{<:AbstractString}`: variable names, one per slice of `sdata`
    along dimension 2.
- `accnames::AbstractVector`: `String`s and/or `Regex`es selecting the variables to sum.
- `addname::AbstractString`: name of the new variable; must not be in `vnames`.
- `sdata::AbstractArray`: data with simulations along dimension 1 and variables along
    dimension 2.

# Keywords
- `attrib=nothing`: vector of attribute dictionaries, one per variable. The new variable
    gets a copy of the attributes of the first selected variable, with `"long_name"` set
    to `addlongname`.
- `addlongname::AbstractString=addname`: `"long_name"` attribute of the new variable.

# Returns
- `Tuple`: `(vnames, sdata, attrib)`, each with the new variable appended.

# Throws
- `ArgumentError`: if an entry of `accnames` selects no variable, `addname` already
    exists, or the sizes of `vnames`, `sdata`, and `attrib` do not agree.

See also [`subtract_variables`](@ref), [`divide_variables`](@ref).

# Examples
```jldoctest
julia> sdata = reshape(Float32[1, 2, 3, 10, 20, 30], 1, 3, 2);  # 1 run, 3 vars, 2 times

julia> vnames, sdata, _ = accumulate_variables(
           ["Clay", "Silt1", "Silt2"], ["Clay", "Silt1"], "PM4Dust", sdata
       );

julia> vnames[end], sdata[1, end, :]
("PM4Dust", Float32[3.0, 30.0])
```
"""
function accumulate_variables(
    vnames::AbstractVector{<:AbstractString},
    accnames::AbstractVector,
    addname::AbstractString,
    sdata::AbstractArray;
    attrib=nothing,
    addlongname::AbstractString=addname,
)
    _check_inputs(vnames, sdata, attrib, (addname,))
    isempty(accnames) && throw(ArgumentError("accnames must not be empty"))

    idxs = Int[]
    for name in accnames
        found = _findvars(vnames, name)
        isempty(found) && throw(ArgumentError("variable $(repr(name)) not found in vnames"))
        append!(idxs, found)
    end
    unique!(idxs)

    accdata = sum(selectdim(sdata, VDIM, idxs); dims=VDIM)
    return _append(vnames, sdata, attrib, [addname], [accdata], first(idxs), [addlongname])
end

"""
    subtract_variables(vnames, diffnames, addname, sdata; kwargs...)

Return `vnames`, `sdata`, and `attrib` with a new variable `addname` appended:
`diffnames[1] - diffnames[2] - … - diffnames[end]`.

Each entry of `diffnames` must name exactly one variable (exact match).

# Arguments
- `vnames::AbstractVector{<:AbstractString}`: variable names, one per slice of `sdata`
    along dimension 2.
- `diffnames::AbstractVector{<:AbstractString}`: minuend followed by one or more
    subtrahends.
- `addname::AbstractString`: name of the new variable; must not be in `vnames`.
- `sdata::AbstractArray`: data with simulations along dimension 1 and variables along
    dimension 2.

# Keywords
- `attrib=nothing`: vector of attribute dictionaries, one per variable. The new variable
    gets a copy of the minuend's attributes, with `"long_name"` set to `addlongname`.
- `addlongname::AbstractString=addname`: `"long_name"` attribute of the new variable.

# Returns
- `Tuple`: `(vnames, sdata, attrib)`, each with the new variable appended.

# Throws
- `ArgumentError`: if `diffnames` has fewer than two entries or an entry is not found,
    `addname` already exists, or the sizes of `vnames`, `sdata`, and `attrib` do not agree.

See also [`accumulate_variables`](@ref), [`divide_variables`](@ref).

# Examples
```jldoctest
julia> sdata = reshape(Float32[5, 2, 1, 50, 20, 10], 1, 3, 2);

julia> vnames, sdata, _ = subtract_variables(["a", "b", "c"], ["a", "b", "c"], "d", sdata);

julia> sdata[1, end, :]
2-element Vector{Float32}:
  2.0
 20.0
```
"""
function subtract_variables(
    vnames::AbstractVector{<:AbstractString},
    diffnames::AbstractVector{<:AbstractString},
    addname::AbstractString,
    sdata::AbstractArray;
    attrib=nothing,
    addlongname::AbstractString=addname,
)
    _check_inputs(vnames, sdata, attrib, (addname,))
    idxs = _findoperands(vnames, diffnames, "diffnames")

    diffdata = copy(selectdim(sdata, VDIM, idxs[1]:idxs[1]))
    for i in idxs[2:end]
        diffdata .-= selectdim(sdata, VDIM, i:i)
    end
    return _append(vnames, sdata, attrib, [addname], [diffdata], idxs[1], [addlongname])
end

"""
    divide_variables(vnames, quotnames, addname, sdata; kwargs...)

Return `vnames`, `sdata`, and `attrib` with a new variable `addname` appended:
`quotnames[1] / quotnames[2] / … / quotnames[end]`.

Each entry of `quotnames` must name exactly one variable (exact match). Division by zero
follows IEEE rules (`Inf` or `NaN`).

# Arguments
- `vnames::AbstractVector{<:AbstractString}`: variable names, one per slice of `sdata`
    along dimension 2.
- `quotnames::AbstractVector{<:AbstractString}`: numerator followed by one or more
    denominators.
- `addname::AbstractString`: name of the new variable; must not be in `vnames`.
- `sdata::AbstractArray`: data with simulations along dimension 1 and variables along
    dimension 2.

# Keywords
- `attrib=nothing`: vector of attribute dictionaries, one per variable. The new variable
    gets a copy of the numerator's attributes with `"long_name"` set to `addlongname`,
    `"scale_factor"` set to `1.0`, and `"add_offset"` (if present) set to `0.0`. Its
    `"units"` attribute is copied unchanged and usually needs to be updated.
- `addlongname::AbstractString=addname`: `"long_name"` attribute of the new variable.

# Returns
- `Tuple`: `(vnames, sdata, attrib)`, each with the new variable appended.

# Throws
- `ArgumentError`: if `quotnames` has fewer than two entries or an entry is not found,
    `addname` already exists, or the sizes of `vnames`, `sdata`, and `attrib` do not agree.

See also [`accumulate_variables`](@ref), [`subtract_variables`](@ref).

# Examples
```jldoctest
julia> sdata = reshape(Float32[6, 3, 60, 20], 1, 2, 2);

julia> vnames, sdata, _ = divide_variables(["a", "b"], ["a", "b"], "ratio", sdata);

julia> sdata[1, end, :]
2-element Vector{Float32}:
 2.0
 3.0
```
"""
function divide_variables(
    vnames::AbstractVector{<:AbstractString},
    quotnames::AbstractVector{<:AbstractString},
    addname::AbstractString,
    sdata::AbstractArray;
    attrib=nothing,
    addlongname::AbstractString=addname,
)
    _check_inputs(vnames, sdata, attrib, (addname,))
    idxs = _findoperands(vnames, quotnames, "quotnames")

    quotdata =
        selectdim(sdata, VDIM, idxs[1]:idxs[1]) ./ selectdim(sdata, VDIM, idxs[2]:idxs[2])
    for i in idxs[3:end]
        quotdata ./= selectdim(sdata, VDIM, i:i)
    end

    vnames, sdata, attrib = _append(
        vnames, sdata, attrib, [addname], [quotdata], idxs[1], [addlongname]
    )
    if attrib !== nothing
        newattr = attrib[end]
        newattr["scale_factor"] = 1.0
        haskey(newattr, "add_offset") && (newattr["add_offset"] = 0.0)
    end
    return vnames, sdata, attrib
end

"""
    accumulate_dust_masses_modele(physvar, vnames, sdata; kwargs...)

Return `vnames`, `sdata`, and `attrib` with cumulative dust mass classes appended, built
from the ModelE size bins `Clay`, `Silt1`, …, `Silt5`. The bin variables are named
`<bin><physvar>`, e.g. `"Clay_emission"` for `physvar = "_emission"`.

The classes are `PM2Dust` (`Clay`), `PM4Dust` (`+ Silt1`), `PM8Dust` (`+ Silt2`),
`PM16Dust` (`+ Silt3`), `PM32Dust` (`+ Silt4`), and `PM64Dust` (`+ Silt5`), named
`<class><physvar>`. Accumulation stops at the first missing bin, so only the classes
complete up to that bin are created; if `Clay` is missing, the inputs are returned
unchanged.

# Arguments
- `physvar::AbstractString`: suffix of the bin variable names, e.g. `"_emission"`.
- `vnames::AbstractVector{<:AbstractString}`: variable names, one per slice of `sdata`
    along dimension 2.
- `sdata::AbstractArray`: data with simulations along dimension 1 and variables along
    dimension 2.

# Keywords
- `attrib=nothing`: vector of attribute dictionaries, one per variable. Each class gets a
    copy of the `Clay` attributes with `"long_name"` set to `"<class> <addlongname>"`.
- `addlongname::AbstractString=""`: text appended to the class name in `"long_name"`.

# Returns
- `Tuple`: `(vnames, sdata, attrib)`, each with the new classes appended.

# Throws
- `ArgumentError`: if a class variable already exists, or the sizes of `vnames`, `sdata`,
    and `attrib` do not agree.

See also [`accumulate_dust_optics_modele`](@ref).

# Examples
```jldoctest
julia> vnames = ["Clay_emission", "Silt1_emission", "Silt2_emission"];

julia> sdata = reshape(Float32[1, 2, 4], 1, 3);

julia> vnames, sdata, _ = accumulate_dust_masses_modele("_emission", vnames, sdata);

julia> vnames[4:end], sdata[1, 4:end]
(["PM2Dust_emission", "PM4Dust_emission", "PM8Dust_emission"], Float32[1.0, 3.0, 7.0])
```
"""
function accumulate_dust_masses_modele(
    physvar::AbstractString,
    vnames::AbstractVector{<:AbstractString},
    sdata::AbstractArray;
    attrib=nothing,
    addlongname::AbstractString="",
)
    srcnames = [bin * physvar for bin in MASS_BINS]
    newnames = [class * physvar for class in MASS_CLASSES]
    return _accumulate_sequence(
        vnames, sdata, attrib, srcnames, newnames, MASS_CLASSES, addlongname
    )
end

"""
    accumulate_dust_optics_modele(physvar, vnames, sdata; kwargs...)

Return `vnames`, `sdata`, and `attrib` with cumulative dust optics classes appended, built
from the ModelE size bins `Clay1`, …, `Clay4`, `Silt1`, …, `Silt5`. The bin variables are
named `<physvar><bin>`, e.g. `"ExtClay1"` for `physvar = "Ext"`.

The classes are `SumClay1Dust` (`Clay1`), `SumClay2Dust` (`+ Clay2`), `SumClay3Dust`
(`+ Clay3`), `PM2Dust` (`+ Clay4`), `PM4Dust` (`+ Silt1`), …, `PM64Dust` (`+ Silt5`), named
`<physvar><class>`. Accumulation stops at the first missing bin, so only the classes
complete up to that bin are created; if `Clay1` is missing, the inputs are returned
unchanged.

# Arguments
- `physvar::AbstractString`: prefix of the bin variable names, e.g. `"Ext"` or `"Abs"`.
- `vnames::AbstractVector{<:AbstractString}`: variable names, one per slice of `sdata`
    along dimension 2.
- `sdata::AbstractArray`: data with simulations along dimension 1 and variables along
    dimension 2.

# Keywords
- `attrib=nothing`: vector of attribute dictionaries, one per variable. Each class gets a
    copy of the `Clay1` attributes with `"long_name"` set to `"<class> <addlongname>"`.
- `addlongname::AbstractString=""`: text appended to the class name in `"long_name"`.

# Returns
- `Tuple`: `(vnames, sdata, attrib)`, each with the new classes appended.

# Throws
- `ArgumentError`: if a class variable already exists, or the sizes of `vnames`, `sdata`,
    and `attrib` do not agree.

See also [`accumulate_dust_masses_modele`](@ref).

# Examples
```jldoctest
julia> vnames = ["ExtClay1", "ExtClay2"];

julia> sdata = reshape(Float32[0.1, 0.2], 1, 2);

julia> vnames, sdata, _ = accumulate_dust_optics_modele("Ext", vnames, sdata);

julia> vnames[3:end]
2-element Vector{String}:
 "ExtSumClay1Dust"
 "ExtSumClay2Dust"
```
"""
function accumulate_dust_optics_modele(
    physvar::AbstractString,
    vnames::AbstractVector{<:AbstractString},
    sdata::AbstractArray;
    attrib=nothing,
    addlongname::AbstractString="",
)
    srcnames = [physvar * bin for bin in OPTICS_BINS]
    newnames = [physvar * class for class in OPTICS_CLASSES]
    return _accumulate_sequence(
        vnames, sdata, attrib, srcnames, newnames, OPTICS_CLASSES, addlongname
    )
end

# Indices of the variables selected by an exact name or a pattern.
_findvars(vnames, name::AbstractString) = findall(==(name), vnames)
_findvars(vnames, pattern::Regex) = findall(v -> occursin(pattern, v), vnames)

# Indices of the operands of a subtraction or division, each found by exact name.
function _findoperands(vnames, names, argname)
    if length(names) < 2
        throw(ArgumentError("$argname needs at least two names, got $(length(names))"))
    end
    idxs = Int[]
    for name in names
        idx = findfirst(==(name), vnames)
        idx === nothing &&
            throw(ArgumentError("variable $(repr(name)) not found in vnames"))
        push!(idxs, idx)
    end
    return idxs
end

# Validate the shared inputs and that none of the new names exists already.
function _check_inputs(vnames, sdata, attrib, newnames)
    ndims(sdata) >= VDIM ||
        throw(ArgumentError("sdata must have at least $VDIM dimensions"))
    if size(sdata, VDIM) != length(vnames)
        nv = size(sdata, VDIM)
        msg = "size(sdata, $VDIM) = $nv but length(vnames) = $(length(vnames))"
        throw(ArgumentError(msg))
    end
    if attrib !== nothing && length(attrib) != length(vnames)
        msg = "length(attrib) = $(length(attrib)) but length(vnames) = $(length(vnames))"
        throw(ArgumentError(msg))
    end
    for name in newnames
        name in vnames && throw(ArgumentError("variable $(repr(name)) already exists"))
    end
    return nothing
end

# Append new variables (names, data slices, attributes copied from variable `srcidx`)
# without modifying the inputs.
function _append(vnames, sdata, attrib, newnames, slices, srcidx, longnames)
    newvnames = String[vnames; newnames]
    newsdata = cat(sdata, slices...; dims=VDIM)
    newattrib = attrib === nothing ? nothing : deepcopy(attrib)
    if newattrib !== nothing
        for longname in longnames
            attr = deepcopy(attrib[srcidx])
            attr["long_name"] = longname
            push!(newattrib, attr)
        end
    end
    return newvnames, newsdata, newattrib
end

# Running sums over `srcnames` in order, one class per bin, stopping at the first missing
# bin. Classes get the attributes of the first bin.
function _accumulate_sequence(vnames, sdata, attrib, srcnames, newnames, classes, longtext)
    _check_inputs(vnames, sdata, attrib, ())
    slices = AbstractArray[]
    names = String[]
    longnames = String[]
    total = nothing
    for (src, new, class) in zip(srcnames, newnames, classes)
        idx = findfirst(==(src), vnames)
        idx === nothing && break
        new in vnames && throw(ArgumentError("variable $(repr(new)) already exists"))
        slice = selectdim(sdata, VDIM, idx:idx)
        total = total === nothing ? copy(slice) : total .+ slice
        push!(slices, total)
        push!(names, new)
        push!(longnames, isempty(longtext) ? class : class * " " * longtext)
    end
    if isempty(slices)
        return String.(vnames), sdata, attrib === nothing ? nothing : deepcopy(attrib)
    end
    srcidx = findfirst(==(srcnames[1]), vnames)
    return _append(vnames, sdata, attrib, names, slices, srcidx, longnames)
end

end # module DerivePhysicsVariablesModelE
