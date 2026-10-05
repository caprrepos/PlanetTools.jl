module Constructors

export fnames_modelE_2d, compose_palette

using ColorSchemes, Dates, Makie

function fnames_modelE_2d(fyear::Int64, fmonth::Int64, lyear::Int64, lmonth::Int64,
        gname::String, ids::AbstractArray{String, 1}; ext::String = ".nc")
# returns a string array with file names of ModelE monthly files over a list of
# run IDs and date ranges
# input:
# fyear:  first year
# fmonth: first month
# lyear:  last year
# lmonth: lmonth
# gname:  file group name (e.g. "aij")
# ids:    run ID list
# ext:    file extension (defaults to ".nc")
# output: 2 dimensional array with file names for (no. run IDs, no. date range)

    nsims = length(ids)
    dvar = Date(fyear, fmonth):Month(1):Date(lyear, lmonth)
    ndates = length(dvar)
    cmy = Array{String}(undef, ndates)

    for i in eachindex(dvar)
        cy = string(dvar[i])[1:4]
        cm = uppercase(monthabbr(dvar[i]))
        cmy[i] = cm*cy
    end

    fnames = Array{String}(undef, nsims, ndates)
    for nr in eachindex(ids)
        fnames[nr,:] = cmy[:] .* "." .* gname .* ids[nr] .* ext
    end
    return fnames
end # fnames_modelE_2d

function compose_palette(breaks; palette1 = nothing, midcolor = nothing,
        palette2 = nothing, demarcation = nothing)
# composes a new palette from one or two input palettes and an optional midcolor.
# If there is only one palette in the input arguments, the midcolor is put as the
# first or last color, depending on whether palette1 or palette2 is given in the
# input arguments

    if length(breaks) <= 1
        println("Error in composePalette: Number of breaks must be greater > 1")
        return
    end

    if demarcation === nothing
        demarcation = 1.
    end

    if palette1 != nothing && palette2 != nothing
        c1 = count(<=(demarcation), breaks) - 1
        c3 = count(>(demarcation), breaks)
        if midcolor != nothing
            c3 -= 1
        end
    elseif palette1 != nothing
        c1 = length(breaks) - 1
        if midcolor != nothing
            c1 -= 1
        end
        c3 = 0
    elseif palette2 != nothing
        c1 = 0
        c3 = length(breaks) - 1
        if midcolor != nothing
            c3 -= 1
        end
    end

    returnPalette = Vector{Any}[]

    if c1 === 1
        lc = Int(length(colorschemes[palette1]) / 2)
        cols1 = (colorschemes[palette1])[lc]
    elseif c1 > 1
        cols1 = ColorScheme(get(colorschemes[palette1], range(0., 1., length = c1))).colors
    end
    if c1 >= 1
        returnPalette = vcat(returnPalette, cols1)
    end

    if midcolor != nothing && c1 >= 0 && c3 >= 0
        cols2 = cgrad([midcolor])[1]
        returnPalette = vcat(returnPalette, cols2)
    end

    if c3 === 1
        lc = Int(length(colorschemes[palette2]) / 2)
        cols3 = (colorschemes[palette2])[lc]
    elseif c3 > 1
        cols3 = ColorScheme(get(colorschemes[palette2], range(0., 1., length = c3))).colors
    end
    if c3 >= 1
        returnPalette = vcat(returnPalette, cols3)
    end

    return returnPalette
end # compose_palette

end # module Constructors
