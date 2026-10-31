## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

# Tests for `DerivePhysicsVariablesModelE`. Imports are in `runtests.jl`.

# Data with value 100 * variable index + linear index of the remaining dimensions, so that
# every variable and grid point is distinct.
function testdata(T, dims...)
    data = Array{T}(undef, dims...)
    for I in CartesianIndices(data)
        data[I] = 100 * I[2] + LinearIndices(data)[I] / length(data)
    end
    return data
end

# Variable slices along dimension 2, keeping that dimension.
vslice(sdata, i::Integer) = selectdim(sdata, 2, i:i)
vslice(sdata, idxs::AbstractVector) = selectdim(sdata, 2, idxs)

attribs(vnames) = [Dict{String,Any}("long_name" => n, "units" => "u_$n") for n in vnames]

@testset "accumulate_variables" begin
    vnames = ["Clay", "Silt1", "Silt2", "SumSilt1"]
    sdata = testdata(Float32, 2, 4, 3, 2)
    attrib = attribs(vnames)

    @testset "sum, shape, and metadata" begin
        v2, s2, a2 = accumulate_variables(
            vnames, ["Clay", "Silt1"], "PM4Dust", sdata; attrib=attrib
        )
        @test v2 == [vnames; "PM4Dust"]
        @test size(s2) == (2, 5, 3, 2)
        @test vslice(s2, 1:4) == sdata
        @test vslice(s2, 5) ≈ vslice(sdata, 1) .+ vslice(sdata, 2)
        @test a2[1:4] == attrib
        @test a2[5] == Dict("long_name" => "PM4Dust", "units" => "u_Clay")
        _, _, a3 = accumulate_variables(
            vnames, ["Silt1", "Clay"], "X", sdata; attrib=attrib, addlongname="Sum"
        )
        @test a3[5] == Dict("long_name" => "Sum", "units" => "u_Silt1")  # First selected.
    end

    @testset "inputs are not modified" begin
        vcopy, scopy, acopy = copy(vnames), copy(sdata), deepcopy(attrib)
        accumulate_variables(vnames, ["Clay"], "X", sdata; attrib=attrib)
        @test vnames == vcopy && sdata == scopy && attrib == acopy
    end

    @testset "any number of dimensions" begin
        for dims in ((2, 4), (2, 4, 3), (2, 4, 3, 2), (1, 4, 2, 2, 3), (1, 4, 2, 2, 3, 2))
            s = testdata(Float64, dims...)
            _, s2, _ = accumulate_variables(vnames, ["Clay", "Silt2"], "X", s)
            @test size(s2) == Base.setindex(dims, 5, 2)
            @test vslice(s2, 5) ≈ vslice(s, 1) .+ vslice(s, 3)
        end
    end

    @testset "name matching" begin
        # Exact names: "Silt1" must not also select "SumSilt1".
        _, s2, _ = accumulate_variables(vnames, ["Silt1"], "X", sdata)
        @test vslice(s2, 5) == vslice(sdata, 2)
        # Regex selects several variables; a variable selected twice is summed once.
        _, s2, _ = accumulate_variables(vnames, [r"^Silt", "Silt1"], "X", sdata)
        @test vslice(s2, 5) ≈ vslice(sdata, 2) .+ vslice(sdata, 3)
        # SubString names.
        v2, _, _ = accumulate_variables(
            split("Clay Silt1 Silt2 SumSilt1"), ["Clay"], "X", sdata
        )
        @test v2 isa Vector{String}
    end

    @testset "element types" begin
        s64 = fill(1 + 1e-9, 1, 2)
        @test accumulate_variables(["a", "b"], ["a", "b"], "c", s64)[2][1, 3] == 2 + 2e-9
        smiss = Union{Missing,Float32}[1 missing]
        @test ismissing(accumulate_variables(["a", "b"], ["a", "b"], "c", smiss)[2][1, 3])
    end

    @testset "invalid input" begin
        @test_throws ArgumentError accumulate_variables(vnames, ["NoSuch"], "X", sdata)
        @test_throws ArgumentError accumulate_variables(
            vnames, ["Clay", "NoSuch"], "X", sdata
        )
        @test_throws ArgumentError accumulate_variables(vnames, [r"^NoSuch"], "X", sdata)
        @test_throws ArgumentError accumulate_variables(vnames, String[], "X", sdata)
        @test_throws ArgumentError accumulate_variables(vnames, ["Clay"], "Silt1", sdata)
        @test_throws ArgumentError accumulate_variables(vnames[1:3], ["Clay"], "X", sdata)
        @test_throws ArgumentError accumulate_variables(
            vnames, ["Clay"], "X", sdata; attrib=attrib[1:3]
        )
        @test_throws ArgumentError accumulate_variables(["a"], ["a"], "X", [1.0])
    end
end

@testset "subtract_variables" begin
    vnames = ["a", "b", "c", "ab"]
    sdata = testdata(Float32, 2, 4, 3)
    attrib = attribs(vnames)

    v2, s2, a2 = subtract_variables(vnames, ["a", "b", "c"], "d", sdata; attrib=attrib)
    @test v2 == [vnames; "d"]
    @test vslice(s2, 5) ≈ vslice(sdata, 1) .- vslice(sdata, 2) .- vslice(sdata, 3)
    @test a2[5] == Dict("long_name" => "d", "units" => "u_a")
    @test vslice(subtract_variables(vnames, ["c", "a"], "d", sdata)[2], 5) ≈
        vslice(sdata, 3) .- vslice(sdata, 1)

    s6 = testdata(Float64, 1, 4, 2, 2, 3, 2)
    @test vslice(subtract_variables(vnames, ["b", "a"], "d", s6)[2], 5) ≈
        vslice(s6, 2) .- vslice(s6, 1)

    vcopy, scopy = copy(vnames), copy(sdata)
    subtract_variables(vnames, ["a", "b"], "d", sdata)
    @test vnames == vcopy && sdata == scopy

    @test_throws ArgumentError subtract_variables(vnames, ["a"], "d", sdata)
    @test_throws ArgumentError subtract_variables(vnames, ["a", "NoSuch"], "d", sdata)
    @test_throws ArgumentError subtract_variables(vnames, ["NoSuch", "a"], "d", sdata)
    @test_throws ArgumentError subtract_variables(vnames, ["a", "b"], "ab", sdata)
end

@testset "divide_variables" begin
    vnames = ["a", "b", "c"]
    sdata = testdata(Float32, 2, 3, 3)
    attrib = [
        Dict{String,Any}("long_name" => n, "scale_factor" => 2.0, "add_offset" => 5.0) for
        n in vnames
    ]

    v2, s2, a2 = divide_variables(vnames, ["a", "b", "c"], "q", sdata; attrib=attrib)
    @test v2 == [vnames; "q"]
    @test vslice(s2, 4) ≈ vslice(sdata, 1) ./ vslice(sdata, 2) ./ vslice(sdata, 3)
    @test a2[4] == Dict("long_name" => "q", "scale_factor" => 1.0, "add_offset" => 0.0)
    @test a2[1:3] == attrib

    _, _, a3 = divide_variables(vnames, ["a", "b"], "q", sdata; attrib=attribs(vnames))
    @test !haskey(a3[4], "add_offset")
    @test a3[4]["scale_factor"] == 1.0

    @test divide_variables(["a", "b"], ["a", "b"], "q", [6 3])[2] == [6 3 2.0]
    @test divide_variables(["a", "b"], ["a", "b"], "q", [1.0 0.0])[2][1, 3] == Inf

    @test_throws ArgumentError divide_variables(vnames, ["a"], "q", sdata)
    @test_throws ArgumentError divide_variables(vnames, ["a", "NoSuch"], "q", sdata)
    @test_throws ArgumentError divide_variables(vnames, ["a", "b"], "c", sdata)
end

@testset "accumulate_dust_masses_modele" begin
    bins = ["Clay", "Silt1", "Silt2", "Silt3", "Silt4", "Silt5"]
    classes = ["PM2Dust", "PM4Dust", "PM8Dust", "PM16Dust", "PM32Dust", "PM64Dust"]
    vnames = [["OtherVar"]; bins .* "_emission"]
    sdata = testdata(Float64, 2, 7, 2, 3)
    attrib = attribs(vnames)

    @testset "all bins" begin
        v2, s2, a2 = accumulate_dust_masses_modele(
            "_emission", vnames, sdata; attrib=attrib, addlongname="Emission"
        )
        @test v2 == [vnames; classes .* "_emission"]
        @test size(s2) == (2, 13, 2, 3)
        for k in 1:6
            @test vslice(s2, 7 + k) ≈ sum(vslice(sdata, 1 + j) for j in 1:k)
        end
        @test [a["long_name"] for a in a2[8:end]] == classes .* " Emission"
        @test all(a["units"] == "u_Clay_emission" for a in a2[8:end])
    end

    @testset "stops at the first missing bin" begin
        keep = [1, 2, 3, 5, 6, 7]  # Without Silt2.
        v2, s2, _ = accumulate_dust_masses_modele(
            "_emission", vnames[keep], vslice(sdata, keep)
        )
        @test v2[7:end] == ["PM2Dust_emission", "PM4Dust_emission"]
        @test vslice(s2, 8) ≈ vslice(sdata, 2) .+ vslice(sdata, 3)
    end

    @testset "no bins found" begin
        v2, s2, a2 = accumulate_dust_masses_modele("_load", vnames, sdata; attrib=attrib)
        @test v2 == vnames && s2 === sdata && a2 == attrib
    end

    @testset "defaults, types, repeated call" begin
        _, _, a2 = accumulate_dust_masses_modele("_emission", vnames, sdata; attrib=attrib)
        @test a2[8]["long_name"] == "PM2Dust"
        smiss = Union{Missing,Float32}[missing 1 2;;;]
        _, s2, _ = accumulate_dust_masses_modele("", ["Clay", "Silt1", "Silt2"], smiss)
        @test all(ismissing, s2[1, 4:6, 1])
        s64 = fill(1 + 1e-9, 1, 2)
        @test accumulate_dust_masses_modele("", ["Clay", "Silt1"], s64)[2][1, 4] == 2 + 2e-9
        v2, s2, _ = accumulate_dust_masses_modele("_emission", vnames, sdata)
        @test_throws ArgumentError accumulate_dust_masses_modele("_emission", v2, s2)
        @test_throws ArgumentError accumulate_dust_masses_modele(
            "_emission", vnames, vslice(sdata, 1:3)
        )
    end
end

@testset "accumulate_dust_optics_modele" begin
    bins = ["Clay1", "Clay2", "Clay3", "Clay4", "Silt1", "Silt2", "Silt3", "Silt4", "Silt5"]
    classes = [
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
    vnames = "Ext" .* bins
    sdata = testdata(Float32, 1, 9, 2, 2, 3, 2)

    v2, s2, a2 = accumulate_dust_optics_modele(
        "Ext", vnames, sdata; attrib=attribs(vnames), addlongname="Extinction"
    )
    @test v2 == [vnames; "Ext" .* classes]
    for k in 1:9
        @test vslice(s2, 9 + k) ≈ sum(vslice(sdata, j) for j in 1:k)
    end
    @test [a["long_name"] for a in a2[10:end]] == classes .* " Extinction"
    @test all(a["units"] == "u_ExtClay1" for a in a2[10:end])

    v2, _, _ = accumulate_dust_optics_modele(
        "Ext", vnames[[1, 2, 4]], vslice(sdata, [1, 2, 4])
    )
    @test v2[4:end] == ["ExtSumClay1Dust", "ExtSumClay2Dust"]
    @test accumulate_dust_optics_modele("Abs", vnames, sdata)[1] == vnames
end
