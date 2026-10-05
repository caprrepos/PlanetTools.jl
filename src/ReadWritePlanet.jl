module ReadWritePlanet

export readnetcdf, writenetcdf

using DataStructures, NCDatasets

# --- Internal Dispatch Functions ---

# Specific case to resolve ambiguity: when both are nothing
function _get_ndatadims(::Nothing, ::Nothing)
    return 0
end

# Case 1: No data provided, but names might exist
function _get_ndatadims(::Nothing, names::Any)
    return 0
end

# Case 2: Data provided but no names
function _get_ndatadims(data::Any, ::Nothing)
    return ndims(data)
end

# Case 3: Multiple variables packed into one array [Var, Lon, Lat, ...]
function _get_ndatadims(data::AbstractArray, names::Vector{String})
    return ndims(data) - 1
end

# Case 4: A single variable array [Lon, Lat, ...]
function _get_ndatadims(data::AbstractArray, names::String)
    return ndims(data)
end

# Case 5: Scalars (Int, Float, etc.)
function _get_ndatadims(data::Number, names::String)
    return 0
end

# --- Main Writing Function ---

function writenetcdf(path;
        mode = "c",
        format::Symbol = :netcdf4,
        data = nothing,
        vtype = nothing,
        names = nothing,
        datadim1 = nothing,
        dim1attr = nothing,
        datadim2 = nothing,
        dim2attr = nothing,
        datadim3 = nothing,
        dim3attr = nothing,
        datadim4 = nothing,
        dim4attr = nothing,
        datadim5 = nothing,
        dim5attr = nothing,
        datadim6 = nothing,
        dim6attr = nothing,
        dimnames = nothing,
        dataattr = nothing,
        datevar = nothing,
        datevarattr = nothing,
        axypvar = nothing,
        axypvarattr = nothing,
        globalattr = nothing
    )

    # 1. Determine number of dimensions
    ndatadims = _get_ndatadims(data, names)

    # 2. Auto-generate generic dimension names if missing
    if dimnames === nothing && ndatadims > 0
        dimnames = "dimension" .* string.(1:ndatadims)
    end

    # 3. Define helper for dimension creation
    function defineDimNc(ds, dimname, datadim; dimattr=nothing)
        if !haskey(ds, dimname)
            # Use Inf for time, else the actual length of the coordinate array
            d_len = (dimname == "time") ? Inf : length(datadim)
            defDim(ds, dimname, d_len)

            if dimattr != nothing
                defVar(ds, dimname, datadim, (dimname,); attrib=dimattr)
            else
                defVar(ds, dimname, datadim, (dimname,))
            end
        end

        return
    end

    NCDataset(path, mode; format=format, attrib=[]) do ds
        # --- A. Define Dimensions ---
        if dimnames != nothing
            for n in eachindex(dimnames)
                if n === 1 && datadim1 != nothing
                    defineDimNc(ds, dimnames[n], datadim1; dimattr=dim1attr)
                elseif n === 2 && datadim2 != nothing
                    defineDimNc(ds, dimnames[n], datadim2; dimattr=dim2attr)
                elseif n === 3 && datadim3 != nothing
                    defineDimNc(ds, dimnames[n], datadim3; dimattr=dim3attr)
                elseif n === 4 && datadim4 != nothing
                    defineDimNc(ds, dimnames[n], datadim4; dimattr=dim4attr)
                elseif n === 5 && datadim5 != nothing
                    defineDimNc(ds, dimnames[n], datadim5; dimattr=dim5attr)
                elseif n === 6 && datadim6 != nothing
                    defineDimNc(ds, dimnames[n], datadim6; dimattr=dim6attr)
                end
            end
        end

        # --- B. Standardize Types ---
        if vtype != nothing
            if axypvar != nothing; axypvar = convert.(vtype, axypvar); end
            if data != nothing; data = convert.(vtype, data); end
        end

        # --- C. Write Main Variables (Scalars & Arrays) ---
        if names != nothing
            # NetCDF requires a Tuple for dims; () signifies a scalar
            vdims = dimnames === nothing ? () : Tuple(dimnames)

            if names isa Vector
                for n in eachindex(names)
                    if !haskey(ds, names[n])
                        attr = (dataattr != nothing) ? dataattr[n] : []
                        # Robust slicing: pulls the n-th variable slice
                        vdata = (ndims(data) > length(vdims)) ? selectdim(data, 1, n) : data
                        defVar(ds, names[n], vdata, vdims; attrib = attr)
                    end
                end
            else
                # Single variable or scalar (e.g., LS1, Ptop)
                if !haskey(ds, names)
                    attr = (dataattr != nothing) ? dataattr : []
                    defVar(ds, names, data, vdims; attrib = attr)
                end
            end
        end

        # --- D. Write Support Variables ---
        if datevar != nothing && !haskey(ds, "Dates")
            attr = (datevarattr === nothing) ? [] : datevarattr
            defVar(ds, "Dates", string.(datevar), ("time",); attrib = attr)
        end

        if axypvar != nothing && !haskey(ds, "axyp")
            axyp_dims = ndims(axypvar) == 2 ? ("lon", "lat") : (ndims(axypvar) == 1 ? ("lat",) : ())
            attr = (axypvarattr === nothing) ? [] : axypvarattr
            defVar(ds, "axyp", axypvar, axyp_dims; attrib = attr)
        end

        # --- E. Attach Global Metadata ---
        if globalattr != nothing
            for (key, value) in globalattr
                ds.attrib[key] = value
            end
        end
    end

    return nothing
end # writenetcdf

"""
readnetcdf(filepath::AbstractString, var_names::Vector{String})

Reads a NetCDF file, extracts metadata and loads data for multiple pre-selected variables
into Julia arrays.

# Arguments
- `filepath::AbstractString`: The full path to the NetCDF file.
- `var_names::Vector{String}`: A list of variable names to process and read.
    (e.g., ["tau_ext_aer", "ssa_aer", "g_aer"])

# Returns
- A NamedTuple containing:
    - `global_attributes::OrderedDict{String, Any}`: Global attributes of the file.
    - `variable_info::OrderedDict{String, Any}`: Metadata for each requested variable.
    - `variable_data::OrderedDict{String, Any}`: The data for each requested variable
      as a Julia Array.
"""
function readnetcdf(filepath::AbstractString, var_names::Vector{String})
    # Use the NCDatasets.jl 'do' block for safe, automatic file closing
    NCDataset(filepath, "r") do ds

        # 1. Read Global Attributes
        global_attributes = OrderedDict(String(k) => v for (k, v) in ds.attrib)

        # Initialize dictionaries to hold variable info and data
        variable_info = OrderedDict{String, Any}()
        variable_data = OrderedDict{String, Any}()

        # 2. Process Pre-Selected Variables
        for var_name in var_names
            if haskey(ds, var_name)
                # Access the variable handler (CFVariable)
                var = ds[var_name]

                # Extract metadata
                info = Dict{String, Any}(
                    "name"          => String(name(var)),
                    "dimensions"    => dimnames(var), # Tuple of dimension names
                    "dimension_sizes" => size(var),   # Tuple of dimension sizes
                    "attributes"    => OrderedDict(String(k) => v for (k, v) in var.attrib)
               )

                # 3. Read variable data into a Julia Array
                # The Array operator reads the entire variable data into memory with the same
                # dimensions as in the file
                data_array = Array(var)

                # Store the results
                variable_info[var_name] = info
                variable_data[var_name] = data_array
            else
                @warn "Variable '$var_name' not found in file: $filepath"
            end
        end

        # Return the results as a NamedTuple
        return (
            global_attributes = global_attributes,
            variable_info = variable_info,
            variable_data = variable_data
       )
    end
end # readnetcdf

end # module ReadWritePlanet
