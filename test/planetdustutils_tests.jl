## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-24

# Tests for `PlanetDustUtils`. Imports are in `runtests.jl`.

# Kok (2011), Eq. 6: dV/dlnD with c_V = 12.62 µm, D̄s = 3.4 µm, σs = 3.0, λ = 12 µm.
function kok_dvdlnd(d; ds=3.4, sigma=3.0, lambda=12.0, cv=12.62)
    return d / cv * (1 + erf(log(d / ds) / (sqrt(2) * log(sigma)))) * exp(-(d / lambda)^3)
end

# Independent check of the integrals: composite Simpson rule in ln D.
function simpson_lnd(f, a, b; n=2000)
    x = range(log(a), log(b); length=n + 1)
    w = [i == 0 || i == n ? 1 : (isodd(i) ? 4 : 2) for i in 0:n]
    return step(x) / 3 * sum(w .* f.(exp.(x)))
end

@testset "calc_fractional_dust_masses_bft" begin
    fractions = calc_fractional_dust_masses_bft
    bounds = [0.2, 2.0, 4.0, 8.0, 16.0, 32.0]

    @testset "bin integrals" begin
        ref = [simpson_lnd(kok_dvdlnd, bounds[i], bounds[i + 1]) for i in 1:5]
        @test fractions(bounds) ≈ ref ./ sum(ref) rtol = 1e-8
        @test fractions(bounds; normalize=false) ≈ ref rtol = 2e-3  # c_V rounded in paper.
        @test sum(fractions(bounds)) ≈ 1
        @test fractions(bounds) isa Vector{Float64}
        @test fractions([1, 2, 4]) ≈ fractions([1.0, 2.0, 4.0])
    end

    @testset "normalization over all sizes" begin
        @test only(fractions([0, Inf]; normalize=false)) ≈ 1
        whole = fractions([0, 1, 5, 20, Inf]; normalize=false)
        @test sum(whole) ≈ 1
        @test whole ≈ fractions([0, 1, 5, 20, Inf])
        # Splitting a bin conserves its fraction.
        split = fractions([2, 3, 4]; normalize=false)
        @test sum(split) ≈ only(fractions([2, 4]; normalize=false))
        # Kok (2011) gives c_V = 12.62 µm for the default parameters.
        cv = 12.62 * simpson_lnd(kok_dvdlnd, 1e-4, 200.0)
        @test cv ≈ 12.62 rtol = 3e-3
    end

    @testset "parameters" begin
        coarse = fractions(bounds; median_diameter=6.0)
        @test coarse[end] > fractions(bounds)[end]
        @test fractions(bounds; crack_length=6.0)[end] < fractions(bounds)[end]
        ref = [
            simpson_lnd(
                d -> kok_dvdlnd(d; ds=2.0, sigma=2.0, lambda=20.0), bounds[i], bounds[i + 1]
            ) for i in 1:5
        ]
        @test fractions(bounds; median_diameter=2.0, geometric_sd=2.0, crack_length=20.0) ≈
            ref ./ sum(ref) rtol = 1e-8
    end

    @testset "midpoint method reproduces the earlier implementation" begin
        legacy = [0.01764, 0.1333, 0.3808, 0.4653, 0.002926]
        @test fractions(bounds; method=:midpoint) ≈ legacy rtol = 1e-3
        # The midpoint rule is far off for wide bins, e.g. by > 50 % in the clay bin.
        exact = fractions(bounds)
        @test abs(fractions(bounds; method=:midpoint)[1] / exact[1] - 1) > 0.5
    end

    @testset "invalid input" begin
        @test_throws ArgumentError fractions([2.0])
        @test_throws ArgumentError fractions([8.0, 4.0, 2.0])
        @test_throws ArgumentError fractions([1.0, 2.0, 2.0])
        @test_throws ArgumentError fractions([-1.0, 2.0])
        @test_throws ArgumentError fractions([0.0, 2.0]; method=:midpoint)
        @test_throws ArgumentError fractions([1.0, Inf]; method=:midpoint)
        @test_throws ArgumentError fractions(bounds; method=:trapezoid)
        @test_throws ArgumentError fractions(bounds; median_diameter=0)
        @test_throws ArgumentError fractions(bounds; geometric_sd=1.0)
        @test_throws ArgumentError fractions(bounds; crack_length=-1)
        @test_throws ArgumentError fractions([1e4, 2e4])  # Distribution is zero there.
    end
end
