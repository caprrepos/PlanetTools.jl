## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

# Tests for `ReadWritePlanet`. Imports are in `runtests.jl`.

const NCDIR = mktempdir()
ncpath(name) = joinpath(NCDIR, name)
ncdims(path) = NCDataset(ds -> collect(keys(ds.dim)), path)
ncvars(path) = NCDataset(ds -> collect(keys(ds)), path)

const LON = [0.0, 90.0, 180.0, 270.0]
const LAT = [-45.0, 0.0, 45.0]
const TIMES = [1.0, 2.0]
const COORDS = (;
    dimnames=["lon", "lat", "time"], datadim1=LON, datadim2=LAT, datadim3=TIMES
)

@testset "writenetcdf and readnetcdf" begin
    data = rand(MersenneTwister(1), Float32, 2, 4, 3, 2)  # [Var, Lon, Lat, Time]

    @testset "several variables" begin
        p = ncpath("multi.nc")
        writenetcdf(
            p;
            data=data,
            names=["u", "v"],
            COORDS...,
            dim1attr=Dict("units" => "degrees_east"),
            dataattr=[Dict("units" => "m s-1"), Dict("units" => "m s-2")],
            globalattr=Dict(:title => "test", "source" => "ModelE"),
        )
        r = readnetcdf(p, ["u", "v", "lon"])
        @test r.variable_data["u"] == data[1, :, :, :]
        @test r.variable_data["v"] == data[2, :, :, :]
        @test r.variable_data["lon"] == LON
        info = r.variable_info["v"]
        @test info["name"] == "v" && info["dimensions"] == ("lon", "lat", "time")
        @test info["dimension_sizes"] == (4, 3, 2)
        @test info["attributes"]["units"] == "m s-2"
        @test r.variable_info["lon"]["attributes"]["units"] == "degrees_east"
        @test r.global_attributes["title"] == "test"
        @test r.global_attributes["source"] == "ModelE"
        @test NCDataset(ds -> collect(keys(ds.dim)) ∩ ["time"], p) == ["time"]
        @test NCDataset(ds -> NCDatasets.unlimited(ds.dim), p) == ["time"]
    end

    @testset "single variable, scalar, generic dimensions" begin
        p = ncpath("single.nc")
        writenetcdf(p; data=data[1, :, :, :], names="u", COORDS...)
        @test readnetcdf(p, ["u"]).variable_data["u"] == data[1, :, :, :]

        p = ncpath("scalar.nc")
        writenetcdf(p; data=3.5, names="Ptop")
        @test readnetcdf(p, ["Ptop"]).variable_data["Ptop"][] == 3.5

        p = ncpath("generic.nc")
        writenetcdf(p; data=data[1, :, :, :], names="u")
        @test ncdims(p) == ["dimension1", "dimension2", "dimension3"]
        @test ncvars(p) == ["u"]  # No coordinate variables without coordinates.
    end

    @testset "SubString names" begin
        p = ncpath("substring.nc")
        writenetcdf(p; data=data, names=split("u v"), COORDS...)
        r = readnetcdf(p, split("u v"))
        @test collect(keys(r.variable_data)) == ["u", "v"]
    end

    @testset "append mode" begin
        p = ncpath("append.nc")
        writenetcdf(p; data=data[1, :, :, :], names="u", COORDS...)
        writenetcdf(p; mode="a", data=data[2, :, :, :], names="v", COORDS...)
        writenetcdf(p; mode="a", data=zero(data[1, :, :, :]), names="u", COORDS...)
        r = readnetcdf(p, ["u", "v"])
        @test r.variable_data["u"] == data[1, :, :, :]  # Not overwritten.
        @test r.variable_data["v"] == data[2, :, :, :]

        # Dimensions exist, coordinate variables do not.
        p = ncpath("append_nocoords.nc")
        writenetcdf(p; data=data[1, :, :, 1], names="u", dimnames=["lon", "lat"])
        writenetcdf(
            p;
            mode="a",
            data=data[2, :, :, 1],
            names="v",
            dimnames=["lon", "lat"],
            datadim1=LON,
            datadim2=LAT,
        )
        @test readnetcdf(p, ["lon"]).variable_data["lon"] == LON
    end

    @testset "dates, cell areas, type conversion" begin
        p = ncpath("support.nc")
        writenetcdf(
            p;
            data=data[1, :, :, :],
            names="u",
            dimnames=["longitude", "latitude", "time"],
            datadim1=LON,
            datadim2=LAT,
            datadim3=TIMES,
            datevar=[Date(2000, 1), Date(2000, 2)],
            axypvar=ones(4, 3),
            axypvarattr=Dict("units" => "m2"),
            vtype=Float64,
        )
        @test ncdims(p) == ["longitude", "latitude", "time"]  # No extra lon/lat.
        r = readnetcdf(p, ["Dates", "axyp", "u"])
        @test r.variable_data["Dates"] == ["2000-01-01", "2000-02-01"]
        @test r.variable_info["axyp"]["dimensions"] == ("longitude", "latitude")
        @test eltype(r.variable_data["u"]) == Float64
        @test eltype(r.variable_data["axyp"]) == Float64

        p = ncpath("support2.nc")  # No "time" dimension; 1-D cell areas.
        writenetcdf(
            p;
            data=data[1, :, :, 1],
            names="u",
            dimnames=["lon", "lat"],
            datevar=[Date(2000, 1), Date(2000, 2), Date(2000, 3)],
            axypvar=ones(3),
        )
        r = readnetcdf(p, ["Dates", "axyp"])
        @test r.variable_info["Dates"]["dimensions"] == ("lat",)
        @test r.variable_info["axyp"]["dimensions"] == ("lat",)
    end

    @testset "fill values, missing, NaN" begin
        p = ncpath("fill.nc")
        x = Union{Missing,Float32}[1 missing; NaN 4]
        writenetcdf(p; data=x, names="x", dataattr=Dict("_FillValue" => -999.0f0))
        y = readnetcdf(p, ["x"]).variable_data["x"]
        @test eltype(y) == Union{Missing,Float32} && ismissing(y[1, 2]) && isnan(y[2, 1])
        z = readnetcdf(p, ["x"]; missing_as_nan=true).variable_data["x"]
        @test eltype(z) == Float32 && isnan(z[1, 2]) && z[2, 2] == 4

        p = ncpath("fill_conv.nc")
        writenetcdf(
            p; data=x, names="x", vtype=Float64, dataattr=Dict("_FillValue" => -999.0)
        )
        @test eltype(readnetcdf(p, ["x"]).variable_data["x"]) == Union{Missing,Float64}

        p = ncpath("int.nc")
        writenetcdf(p; data=[1 2; 3 4], names="n")
        @test readnetcdf(p, ["n"]; missing_as_nan=true).variable_data["n"] == [1 2; 3 4]
    end

    @testset "missing variable" begin
        p = ncpath("multi.nc")
        r = @test_logs (:warn, r"not found") readnetcdf(p, ["u", "nosuch"])
        @test collect(keys(r.variable_data)) == ["u"]
    end

    @testset "invalid input is rejected before writing" begin
        p = ncpath("invalid.nc")
        @test_throws ArgumentError writenetcdf(p; names="u")
        @test_throws ArgumentError writenetcdf(
            p; data=data, names=["u", "v", "w"], COORDS...
        )
        @test_throws ArgumentError writenetcdf(
            p; data=data, names=["u", "v"], COORDS..., dataattr=[Dict("units" => "m")]
        )
        @test_throws ArgumentError writenetcdf(
            p; data=data[1, :, :, :], names="u", dimnames=["lon", "lat"]
        )
        @test_throws ArgumentError writenetcdf(
            p; data=data[1, :, 1:2, :], names="u", COORDS...
        )
        @test_throws ArgumentError writenetcdf(p; data=3.0, names=["a", "b"])
        @test !isfile(p)
    end
end
