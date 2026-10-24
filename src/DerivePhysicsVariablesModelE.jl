module DerivePhysicsVariablesModelE

# Summary: Functions for calculating derived physics variables from other variables
# Author: Dr Jan P Perlwitz; Climate, Aerosol, and Pollution Research, LLC
# Contact: jan.p.perlwitz@caprllc.com
export accumulate_variables,
    accumulate_dust_masses_modelE,
    accumulate_dust_optics_modelE,
    divide_variables,
    subtract_variables

function accumulate_variables(
    vnames::Vector{String},
    accnames::Vector{String},
    addname::String,
    sdata;
    attrib=nothing,
    addlongname="not set",
)
    # Summary: accumulates data from multiple variables, whose names match the elements of a subset
    #          of variable names.
    # Input: vnames - string vector with variable names.
    #        accnames - string vector with subset of variable names to be accumulated.
    #        addname - string with new name for accumulated variable.
    #        sdata - array with data (at least 2 dimensions, one for number of the
    #		     simulations and the other one for the number of variables).
    #        attrib - array of dictionaries with attributes.
    #	     addlongname - new long physics name for accumulated variables for NetCDF.
    # Output: vnames - variable names with added names after accumulation.
    #         sdata - simulation data with added data after accumulation.
    #         attrib - array with dictionaries with added attributes for new variables.

    ndimssdata = ndims(sdata)
    nsizesdata = collect(size(sdata))
    nvdim = 2
    for nind in eachindex(nsizesdata)
        if nind === nvdim
            nsizesdata[nind] = 1
        else
            nsizesdata[nind] = nsizesdata[nind]
        end
    end

    accdata = Array{Float32}(undef, Tuple(nsizesdata))

    nvFirst = findfirst(==(1), occursin.(accnames[1], vnames)) # this line is useless, as
    # nvFirst is set to nothing below?

    findacc = Int[]
    nvFirst = nothing
    for name in accnames
        cond = occursin.(name, vnames)
        if any(cond)
            findacc = vcat(findacc, findall(x -> x == 1, cond))
            nvFirst = findfirst(==(1), occursin.(name, vnames))
        end
    end

    if nvFirst != nothing
        accdata = sum(sdata[:, findacc, :, :, :]; dims=nvdim)

        push!(vnames, addname)
        sdata = cat(sdata, accdata; dims=nvdim)

        if ! isnothing(attrib)
            newattributes = deepcopy(attrib[nvFirst])
            newattributes["long_name"] = addlongname
            push!(attrib, newattributes)
        end
    end

    return vnames, sdata, attrib
end # accumulate_variables

function subtract_variables(
    vnames::Vector{String},
    diffnames::Vector{String},
    addname::String,
    sdata;
    attrib=nothing,
    addlongname="not set",
)
    # Summary: subtract subsequent multiple variables, whose names match the elements of a subset
    #          of variable names, from the first variable in the subset.
    # Input: vnames - string vector with variable names.
    #        diffnames - string vector with subset of variable names to be subtracted from
    #            each other.
    #        addname - string with new name for variable after subtraction.
    #        sdata - array with data (at least 2 dimensions, one for number of the
    #		     simulations and the other one for the number of variables).
    #        attrib - array of dictionaries with attributes.
    #	     addlongname - new long physics name of variable after subtraction for NetCDF write.
    # Output: vnames - variable names with added names after subtraction.
    #         sdata - simulation data with added data after subtraction.
    #         attrib - array with dictionaries with added attributes for new variables.

    ndimssdata = ndims(sdata)
    nsizesdata = collect(size(sdata))
    nvdim = 2
    for nind in eachindex(nsizesdata)
        if nind === nvdim
            nsizesdata[nind] = 1
        else
            nsizesdata[nind] = nsizesdata[nind]
        end
    end

    diffdata = Array{Float32}(undef, Tuple(nsizesdata))

    findMinuend = nothing
    cond = occursin.(diffnames[1], vnames)
    if any(cond)
        findMinuend = findall(==(1), cond)
    else
        println("Error: minuend ", diffnames[1], " not found")
        return vnames, sdata, attrib
    end

    findSubtrahend = nothing
    cond = occursin.(diffnames[2], vnames)
    if any(cond)
        findSubtrahend = findall(==(1), cond)
        diffdata = sdata[:, findMinuend, :, :, :] .- sdata[:, findSubtrahend, :, :, :]
    else
        println("Error: first subtrahend ", diffnames[2], " not found")
        return vnames, sdata, attrib
    end

    if isless(2, length(diffnames))
        for name in diffnames[3:end]
            findSubtrahend = nothing
            cond = occursin.(name, vnames)
            if any(cond)
                findSubtrahend = findall(==(1), cond)
                diffdata = diffdata .- sdata[:, findSubtrahend, :, :, :]
            else
                println("Error: subsequent subtrahend ", name, " not found")
                return vnames, sdata, attrib
            end
        end
    end

    push!(vnames, addname)
    sdata = cat(sdata, diffdata; dims=nvdim)

    if ! isnothing(attrib)
        newattributes = deepcopy(attrib[findMinuend[1]])
        newattributes["long_name"] = addlongname
        push!(attrib, newattributes)
    end

    return vnames, sdata, attrib
end # subtract_variables

function divide_variables(
    vnames::Vector{String},
    quotnames::Vector{String},
    addname::String,
    sdata;
    attrib=nothing,
    addlongname="not set",
)
    # Summary: calculate the serial quotient of multiple variables that are in a subset with the
    #          first one being the numerator
    # Input: vnames - string vector with variable names.
    #        quotnames - string vector with subset of variable names from which the serial quotient
    #            is calculated
    #        addname - string with new name for variable after division
    #        sdata - array with data (at least 2 dimensions, one for number of the
    #		     simulations and the other one for the number of variables).
    #        attrib - array of dictionaries with attributes.
    #	     addlongname - new long physics name of variable after division for NetCDF write.
    # Output: vnames - variable names with added names after division.
    #         sdata - simulation data with added data after division.
    #         attrib - array with dictionaries with added attributes for new variables.

    ndimssdata = ndims(sdata)
    nsizesdata = collect(size(sdata))
    nvdim = 2
    for nind in eachindex(nsizesdata)
        if nind === nvdim
            nsizesdata[nind] = 1
        else
            nsizesdata[nind] = nsizesdata[nind]
        end
    end

    quotdata = Array{Float32}(undef, Tuple(nsizesdata))

    findNumerator = nothing
    cond = occursin.(quotnames[1], vnames)
    if any(cond)
        findNumerator = findall(==(1), cond)
    else
        println("Error: numerator ", quotnames[1], " not found")
        return vnames, sdata, attrib
    end

    findDenominator = nothing
    cond = occursin.(quotnames[2], vnames)
    if any(cond)
        findDenominator = findall(==(1), cond)
        quotdata = sdata[:, findNumerator, :, :, :] ./ sdata[:, findDenominator, :, :, :]
    else
        println("Error: first denominator ", quotnames[2], " not found")
        return vnames, sdata, attrib
    end

    if isless(2, length(quotnames))
        for name in quotnames[3:end]
            findDenominator = nothing
            cond = occursin.(name, vnames)
            if any(cond)
                findDenominator = findall(==(1), cond)
                quotdata = quotdata ./ sdata[:, findDenominator, :, :, :]
            else
                println("Error: subsequent denominator ", name, " not found")
                return vnames, sdata, attrib
            end
        end
    end

    push!(vnames, addname)
    sdata = cat(sdata, quotdata; dims=nvdim)

    if ! isnothing(attrib)
        newattributes = deepcopy(attrib[findNumerator[1]])
        newattributes["scale_factor"] = 1.0
        newattributes["long_name"] = addlongname
        push!(attrib, newattributes)
    end

    return vnames, sdata, attrib
end # divide_variables

"""
    accumulate_dust_masses_modelE(physvar::String, vnames::Vector{String}, sdata;
                                 attrib = nothing, addlongname = "not set")

Accumulates dust mass variables from discrete ModelE size bins into cumulative
categories (e.g., SumClay, PM2, PM4, etc.). The function is designed to be rank-agnostic,
supporting 5-D (spatial-temporal) or 6-D (spatial-temporal-vertical) datasets.

# Arguments
- `physvar::String`: The appendix of the physics variable to search for (e.g., "_emission").
- `vnames::Vector{String}`: A list of the current variable names present in the `sdata` array.
  This vector is modified in-place to include the names of the new accumulated variables.
- `sdata::AbstractArray`: The multi-dimensional simulation data. The 2nd dimension
  (index 2) must correspond to the variables. Supports any number of dimensions (e.g., 5-D or 6-D).
- `attrib::Union{Nothing, Vector{OrderedDict{String, Any}}}`: A vector of dictionaries
  containing NetCDF attributes. If provided, it is updated with attributes for the new variables.
- `addlongname::String`: A descriptive suffix added to the "long_name" attribute of
  accumulated variables (e.g., "Optical Depth").

# Returns
- `vnames`: The updated vector of variable names including the new cumulative categories.
- `sdata`: The updated data array with the new cumulative variables appended along the 2nd dimension.
- `attrib`: The updated attribute vector (or `nothing`).

# Accumulation Logic
The function searches for the following bins in sequence:
`Clay1 -> Clay2 -> Clay3 -> Clay4 -> Silt1 -> Silt2 -> Silt3 -> Silt4 -> Silt5`.
If a bin in the sequence is missing, the accumulation **stops immediately**, and no
further cumulative variables are created.

# Accumulated Categories Created
1. **SumClay1Dust**: Clay1
2. **SumClay2Dust**: Clay1 + Clay2
3. **SumClay3Dust**: Clay1 + Clay2 + Clay3
4. **PM2Dust**: SumClay3 + Clay4 (Total Clay)
5. **PM4Dust**: PM2 + Silt1
6. **PM8Dust**: PM4 + Silt2
7. **PM16Dust**: PM8 + Silt3
8. **PM32Dust**: PM16 + Silt4
9. **PM64Dust**: PM32 + Silt5
"""
function accumulate_dust_masses_modelE(
    physvar::String, vnames::Vector{String}, sdata; attrib=nothing, addlongname="not set"
)
    nvdim = 2 # The dimension index for variables

    # Define the dust bins to look for and the sequence of accumulation names
    bins = ["Clay", "Silt1", "Silt2", "Silt3", "Silt4", "Silt5"]
    new_names = ["PM2Dust", "PM4Dust", "PM8Dust", "PM16Dust", "PM32Dust", "PM64Dust"]

    # Setup storage for new data slices to avoid repetitive cat() allocations
    new_slices = []

    # Initialize running total based on the spatial/temporal dimensions of sdata
    # This works for any rank (5D, 6D, etc.)
    spatial_size = collect(size(sdata))
    spatial_size[nvdim] = 1
    running_total = zeros(Float32, Tuple(spatial_size))

    nvFirst = 0

    # Process the bins in sequence
    for i in eachindex(bins)
        target_vname = bins[i] * physvar

        # Check if the bin exists in vnames
        idx = findfirst(==(target_vname), vnames)

        if !isnothing(idx)
            if nvFirst == 0
                nvFirst = idx
            end

            # Extract the bin data and add in-place to the running total
            # selectdim makes this rank-agnostic
            bin_data = sum(selectdim(sdata, nvdim, idx:idx); dims=nvdim)
            running_total .+= bin_data

            # Save a snapshot of the accumulation at this stage
            push!(new_slices, copy(running_total))

            # Metadata Update
            push!(vnames, new_names[i] * physvar)
            if !isnothing(attrib)
                new_attr = deepcopy(attrib[nvFirst])
                new_attr["long_name"] = new_names[i] * " " * addlongname
                push!(attrib, new_attr)
            end
        else
            # STOP requirement: If a bin in the sequence is missing,
            # we stop accumulating further categories immediately.
            break
        end
    end

    # Perform a single concatenation at the end for efficiency
    if !isempty(new_slices)
        sdata = cat(sdata, new_slices...; dims=nvdim)
    end

    return vnames, sdata, attrib
end # accumulate_dust_masses_modelE

"""
    accumulate_dust_optics_modelE(physvar::String, vnames::Vector{String}, sdata;
                                 attrib = nothing, addlongname = "not set")

Accumulates dust optical properties from discrete ModelE size bins into cumulative
categories (e.g., SumClay, PM2, PM4, etc.). The function is designed to be rank-agnostic,
supporting 5-D (spatial-temporal) or 6-D (spatial-temporal-vertical) datasets.

# Arguments
- `physvar::String`: The prefix of the physics variable to search for (e.g., "Ext", "Abs").
- `vnames::Vector{String}`: A list of the current variable names present in the `sdata` array.
  This vector is modified in-place to include the names of the new accumulated variables.
- `sdata::AbstractArray`: The multi-dimensional simulation data. The 2nd dimension
  (index 2) must correspond to the variables. Supports any number of dimensions (e.g., 5-D or 6-D).
- `attrib::Union{Nothing, Vector{OrderedDict{String, Any}}}`: A vector of dictionaries
  containing NetCDF attributes. If provided, it is updated with attributes for the new variables.
- `addlongname::String`: A descriptive suffix added to the "long_name" attribute of
  accumulated variables (e.g., "Optical Depth").

# Returns
- `vnames`: The updated vector of variable names including the new cumulative categories.
- `sdata`: The updated data array with the new cumulative variables appended along the 2nd dimension.
- `attrib`: The updated attribute vector (or `nothing`).

# Accumulation Logic
The function searches for the following bins in sequence:
`Clay1 -> Clay2 -> Clay3 -> Clay4 -> Silt1 -> Silt2 -> Silt3 -> Silt4 -> Silt5`.
If a bin in the sequence is missing, the accumulation **stops immediately**, and no
further cumulative variables are created.

# Accumulated Categories Created
1. **SumClay1Dust**: Clay1
2. **SumClay2Dust**: Clay1 + Clay2
3. **SumClay3Dust**: Clay1 + Clay2 + Clay3
4. **PM2Dust**: SumClay3 + Clay4 (Total Clay)
5. **PM4Dust**: PM2 + Silt1
6. **PM8Dust**: PM4 + Silt2
7. **PM16Dust**: PM8 + Silt3
8. **PM32Dust**: PM16 + Silt4
9. **PM64Dust**: PM32 + Silt5
"""
function accumulate_dust_optics_modelE(
    physvar::String, vnames::Vector{String}, sdata; attrib=nothing, addlongname="not set"
)
    nvdim = 2 # The dimension index for variables

    # Define the dust bins to look for and the sequence of accumulation names
    bins = ["Clay1", "Clay2", "Clay3", "Clay4", "Silt1", "Silt2", "Silt3", "Silt4", "Silt5"]
    new_names = [
        "SumClay1Dust",
        "SumClay2Dust",
        "SumClay3Dust",
        "PM2Dust",
        "PM4Dust",
        "PM8Dust",
        "PM16Dust",
        "PM32Dust",
        "PM64Dust",
    ]

    # Setup storage for new data slices to avoid repetitive cat() allocations
    new_slices = []

    # Initialize running total based on the spatial/temporal dimensions of sdata
    # This works for any rank (5D, 6D, etc.)
    spatial_size = collect(size(sdata))
    spatial_size[nvdim] = 1
    running_total = zeros(Float32, Tuple(spatial_size))

    nvFirst = 0

    # Process the bins in sequence
    for i in eachindex(bins)
        target_vname = physvar * bins[i]

        # Check if the bin exists in vnames
        idx = findfirst(==(target_vname), vnames)

        if !isnothing(idx)
            if nvFirst == 0
                nvFirst = idx
            end

            # Extract the bin data and add in-place to the running total
            # selectdim makes this rank-agnostic
            bin_data = sum(selectdim(sdata, nvdim, idx:idx); dims=nvdim)
            running_total .+= bin_data

            # Save a snapshot of the accumulation at this stage
            push!(new_slices, copy(running_total))

            # Metadata Update
            push!(vnames, physvar * new_names[i])
            if !isnothing(attrib)
                new_attr = deepcopy(attrib[nvFirst])
                new_attr["long_name"] = new_names[i] * " " * addlongname
                push!(attrib, new_attr)
            end
        else
            # STOP requirement: If a bin in the sequence is missing,
            # we stop accumulating further categories immediately.
            break
        end
    end

    # Perform a single concatenation at the end for efficiency
    if !isempty(new_slices)
        sdata = cat(sdata, new_slices...; dims=nvdim)
    end

    return vnames, sdata, attrib
end # function accumulate_dust_optics_modelE

end # module DerivePhysicsVariablesModelE
