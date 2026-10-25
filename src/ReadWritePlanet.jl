## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

"""
    PlanetTools.ReadWritePlanet

Submodule of PlanetTools for reading and writing NetCDF files of planetary model output
with NCDatasets.jl: [`writenetcdf`](@ref) writes variables with coordinates and metadata,
[`readnetcdf`](@ref) reads selected variables and their metadata. Load it with
`using PlanetTools.ReadWritePlanet`.
"""
module ReadWritePlanet

using NCDatasets: NCDataset, defDim, defVar, dimnames, name
using OrderedCollections: OrderedDict

export readnetcdf, writenetcdf

# Maximum number of dimensions with coordinate keywords `datadimN` / `dimNattr`.
const MAX_COORDS = 6

"""
    writenetcdf(path; kwargs...) -> Nothing

Write one or several variables, their coordinates, and metadata to the NetCDF file `path`.

Dimension `i` is named `dimnames[i]` (default `"dimension<i>"`). If a coordinate
`datadim<i>` is given, the dimension gets its length and a coordinate variable of the same
name; a dimension named `"time"` is unlimited. Dimensions without coordinates are created
from the size of the data. Existing dimensions and variables (e.g. with `mode="a"`) are
kept; a variable that already exists is not overwritten.

# Keywords
- `data=nothing`: a scalar, an array for one variable, or for several variables an array
    whose first dimension indexes `names`, `[Var, Lon, Lat, …]`.
- `names=nothing`: variable name (`AbstractString`) or vector of names.
- `dimnames=nothing`: dimension names of the variables (without the `Var` dimension).
- `datadim1`, …, `datadim6`: coordinates of dimensions 1 to 6; `dim1attr`, …, `dim6attr`:
    their attributes (e.g. `Dict("units" => "degrees_east")`).
- `dataattr=nothing`: attributes of the variable, or a vector with one dictionary per name.
- `vtype=nothing`: element type to convert `data` and `axypvar` to, e.g. `Float32`.
- `datevar=nothing`: dates written as strings to variable `"Dates"` along the `"time"`
    dimension (or the last dimension); `datevarattr`: its attributes.
- `axypvar=nothing`: cell areas `(Lon, Lat)` (or `(Lat,)`) written to variable `"axyp"` on
    the first two dimensions (default `"lon"`, `"lat"`); `axypvarattr`: its attributes.
- `globalattr=nothing`: global attributes, any collection of name–value pairs.
- `mode::AbstractString="c"`: `"c"` creates (overwrites) the file, `"a"` appends.
- `format::Symbol=:netcdf4`: file format, see `NCDatasets.NCDataset`.

# Throws
- `ArgumentError`: if the sizes of `data`, `names`, `dimnames`, the coordinates, or
    `dataattr` do not agree. The inputs are checked before the file is opened.

See also [`readnetcdf`](@ref).

# Examples
```jldoctest
julia> path = joinpath(mktempdir(), "example.nc");

julia> u = reshape(Float32.(1:6), 3, 2);  # (lon, lat)

julia> writenetcdf(path; data=u, names="u", dimnames=["lon", "lat"],
           datadim1=[0.0, 120.0, 240.0], datadim2=[-45.0, 45.0],
           dataattr=Dict("units" => "m s-1"), globalattr=Dict("title" => "test"));

julia> r = readnetcdf(path, ["u"]);

julia> r.variable_data["u"] == u, r.variable_info["u"]["dimensions"]
(true, ("lon", "lat"))
```
"""
function writenetcdf(
    path::AbstractString;
    mode::AbstractString="c",
    format::Symbol=:netcdf4,
    data=nothing,
    vtype=nothing,
    names=nothing,
    datadim1=nothing,
    dim1attr=nothing,
    datadim2=nothing,
    dim2attr=nothing,
    datadim3=nothing,
    dim3attr=nothing,
    datadim4=nothing,
    dim4attr=nothing,
    datadim5=nothing,
    dim5attr=nothing,
    datadim6=nothing,
    dim6attr=nothing,
    dimnames=nothing,
    dataattr=nothing,
    datevar=nothing,
    datevarattr=nothing,
    axypvar=nothing,
    axypvarattr=nothing,
    globalattr=nothing,
)
    coords = (datadim1, datadim2, datadim3, datadim4, datadim5, datadim6)
    coordattrs = (dim1attr, dim2attr, dim3attr, dim4attr, dim5attr, dim6attr)
    multi = names isa AbstractVector
    nvdims = _nvardims(data, multi)
    if dimnames === nothing && nvdims > 0
        dimnames = ["dimension$i" for i in 1:nvdims]
    end
    _check_write_inputs(data, names, multi, nvdims, dimnames, coords, dataattr)

    if vtype !== nothing
        data = data === nothing ? nothing : _convert(vtype, data)
        axypvar = axypvar === nothing ? nothing : _convert(vtype, axypvar)
    end

    NCDataset(path, mode; format=format) do ds
        if dimnames !== nothing
            for (i, dname) in enumerate(dimnames)
                i <= MAX_COORDS && coords[i] !== nothing || continue
                _define_coordinate!(ds, dname, coords[i], coordattrs[i])
            end
        end

        if names !== nothing
            vdims = dimnames === nothing ? () : Tuple(dimnames)
            for (n, vname) in enumerate(multi ? names : (names,))
                haskey(ds, vname) && continue
                vdata = multi ? selectdim(data, 1, n) : data
                attr = if multi
                    _attr(dataattr === nothing ? nothing : dataattr[n])
                else
                    _attr(dataattr)
                end
                defVar(ds, vname, vdata, vdims; attrib=attr)
            end
        end

        if datevar !== nothing && !haskey(ds, "Dates")
            tdim = haskey(ds.dim, "time") || dimnames === nothing ? "time" : last(dimnames)
            defVar(ds, "Dates", string.(datevar), (tdim,); attrib=_attr(datevarattr))
        end

        if axypvar !== nothing && !haskey(ds, "axyp")
            horiz = if dimnames === nothing || length(dimnames) < 2
                ("lon", "lat")
            else
                (dimnames[1], dimnames[2])
            end
            adims = if ndims(axypvar) == 2
                horiz
            elseif ndims(axypvar) == 1
                (horiz[2],)
            else
                ()
            end
            defVar(ds, "axyp", axypvar, adims; attrib=_attr(axypvarattr))
        end

        if globalattr !== nothing
            for (key, value) in globalattr
                ds.attrib[String(key)] = value
            end
        end
    end
    return nothing
end

"""
    readnetcdf(filepath, var_names; missing_as_nan=false) -> NamedTuple

Read the global attributes and the selected variables of a NetCDF file, with their metadata.
Values are read as NCDatasets returns them: scale factors and offsets are applied, and
values equal to `_FillValue` become `missing`. Requested variables that do not exist are
skipped with a warning.

# Arguments
- `filepath::AbstractString`: path of the NetCDF file.
- `var_names::AbstractVector{<:AbstractString}`: names of the variables to read.

# Keywords
- `missing_as_nan::Bool=false`: return floating-point data with `NaN` instead of
    `missing` (element type without `Missing`), as expected by PlanetTools' statistics.

# Returns
- `NamedTuple` with fields
    - `global_attributes::OrderedDict{String,Any}`: global attributes;
    - `variable_info::OrderedDict{String,Any}`: per variable, a dictionary with `"name"`,
        `"dimensions"` (tuple of names), `"dimension_sizes"`, and `"attributes"`;
    - `variable_data::OrderedDict{String,Any}`: per variable, its data as an `Array`.

See also [`writenetcdf`](@ref).

# Examples
```jldoctest
julia> path = joinpath(mktempdir(), "example.nc");

julia> x = [1.0 NaN; 3.0 4.0];

julia> writenetcdf(path; data=x, names="x", dataattr=Dict("_FillValue" => -999.0));

julia> eltype(readnetcdf(path, ["x"]).variable_data["x"])
Union{Missing, Float64}

julia> isnan(readnetcdf(path, ["x"]; missing_as_nan=true).variable_data["x"][1, 2])
true
```
"""
function readnetcdf(
    filepath::AbstractString,
    var_names::AbstractVector{<:AbstractString};
    missing_as_nan::Bool=false,
)
    return NCDataset(filepath, "r") do ds
        global_attributes = OrderedDict{String,Any}(String(k) => v for (k, v) in ds.attrib)
        variable_info = OrderedDict{String,Any}()
        variable_data = OrderedDict{String,Any}()
        for vname in var_names
            if !haskey(ds, vname)
                @warn "Variable '$vname' not found in file: $filepath"
                continue
            end
            var = ds[vname]
            variable_info[vname] = OrderedDict{String,Any}(
                "name" => String(name(var)),
                "dimensions" => dimnames(var),
                "dimension_sizes" => size(var),
                "attributes" =>
                    OrderedDict{String,Any}(String(k) => v for (k, v) in var.attrib),
            )
            x = Array(var)
            variable_data[vname] = missing_as_nan ? _missing_to_nan(x) : x
        end
        return (; global_attributes, variable_info, variable_data)
    end
end

# --- Internal helpers ---

# Number of dimensions of each variable in `data`.
_nvardims(::Nothing, multi) = 0
_nvardims(::Number, multi) = 0
_nvardims(data::AbstractArray, multi) = multi ? ndims(data) - 1 : ndims(data)

function _check_write_inputs(data, names, multi, nvdims, dimnames, coords, dataattr)
    if names !== nothing && data === nothing
        throw(ArgumentError("names given without data"))
    end
    if multi
        data isa AbstractArray || throw(ArgumentError("several names need an array"))
        if size(data, 1) != length(names)
            msg = "size(data, 1) = $(size(data, 1)) but there are $(length(names)) names"
            throw(ArgumentError(msg))
        end
        if dataattr !== nothing && length(dataattr) != length(names)
            throw(ArgumentError("dataattr must have one entry per name"))
        end
    end
    if dimnames !== nothing && names !== nothing && length(dimnames) != nvdims
        msg = "$(length(dimnames)) dimnames for variables with $nvdims dimensions"
        throw(ArgumentError(msg))
    end
    n = dimnames === nothing ? 0 : length(dimnames)
    for i in 1:min(n, MAX_COORDS)
        c = coords[i]
        if c !== nothing && data isa AbstractArray && names !== nothing
            len = size(data, multi ? i + 1 : i)
            if length(c) != len
                msg = "datadim$i has length $(length(c)), data dimension $i has size $len"
                throw(ArgumentError(msg))
            end
        end
    end
    return nothing
end

# Define dimension `dname` and its coordinate variable, unless they exist.
function _define_coordinate!(ds, dname, coord, attr)
    if !haskey(ds.dim, dname)
        defDim(ds, dname, dname == "time" ? Inf : length(coord))
    end
    if !haskey(ds, dname)
        defVar(ds, dname, coord, (dname,); attrib=_attr(attr))
    end
    return nothing
end

_attr(::Nothing) = Pair{String,Any}[]
_attr(a) = a

_convert(::Type{T}, x::Number) where {T} = convert(T, x)
function _convert(::Type{T}, x::AbstractArray) where {T}
    return map(v -> v === missing ? missing : convert(T, v), x)
end

_missing_to_nan(x) = x
function _missing_to_nan(x::AbstractArray{Union{Missing,T}}) where {T<:AbstractFloat}
    return T[ismissing(v) ? T(NaN) : v for v in x]
end

end # module ReadWritePlanet
