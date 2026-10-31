## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

# Tests for `PlanetAerosolUtils`. Imports are in `runtests.jl`.

@testset "calc_size_integrated_effective_radius" begin
    reff = calc_size_integrated_effective_radius

    @testset "agrees with the moment definition" begin
        # Monodisperse bins: N particles of radius r and density rho per bin.
        r = [0.1, 0.5, 2.0, 8.0]
        N = [1e6, 1e3, 1.0, 1e-2]
        rho = [2650.0, 2650.0, 1500.0, 3000.0]
        mass = N .* (4 / 3) .* π .* r .^ 3 .* rho
        @test reff(mass, rho, r) ≈ sum(N .* r .^ 3) / sum(N .* r .^ 2)
    end

    @testset "units and special cases" begin
        r = [0.2, 1.0, 4.0]
        mass = [1.0, 3.0, 2.0]  # kg
        rho = 2650.0  # kg m⁻³
        # Same result with mass in g and density in g cm⁻³.
        @test reff(1e3 .* mass, rho / 1e3, r) ≈ reff(mass, rho, r)
        @test reff(mass, rho, r) == reff(mass, fill(rho, 3), r)
        @test reff([5.0], 2650.0, [1.5]) ≈ 1.5
        @test reff(mass, rho, 1e-6 .* r) ≈ 1e-6 * reff(mass, rho, r)  # µm to m.
        @test isnan(reff(zeros(3), rho, r))
        @test reff(Float32.(mass), 2650.0f0, Float32.(r)) isa Float32
    end

    @testset "invalid input" begin
        @test_throws ArgumentError reff(Float64[], 2650.0, Float64[])
        @test_throws ArgumentError reff([1.0, 2.0], 2650.0, [1.0])
        @test_throws ArgumentError reff([1.0, 2.0], [2650.0], [1.0, 2.0])
        @test_throws ArgumentError reff([1.0, 2.0], 2650.0, [0.0, 1.0])
        @test_throws ArgumentError reff([1.0, 2.0], [2650.0, -1.0], [1.0, 2.0])
    end
end

@testset "interpolate_extinction_coeff_log" begin
    interp = interpolate_extinction_coeff_log
    bounds = [(0.2, 0.4), (0.4, 0.8), (0.8, 1.6), (1.6, 3.2)]
    centers = [sqrt(b[1] * b[2]) for b in bounds]
    radii = [0.5, 1.0, 2.0]
    # A power law is linear in log-log space, so interpolation and linear extrapolation
    # must reproduce it exactly.
    powerlaw(λ, r) = 3.0 * λ^-1.5 * r^-0.8
    table = [powerlaw(λ, r) for λ in centers, r in radii]

    @testset "values" begin
        @test interp(bounds, radii, table, centers, radii) ≈ table
        tl, tr = [0.33, 0.7, 2.5], [0.6, 1.7]
        @test interp(bounds, radii, table, tl, tr) ≈ [powerlaw(λ, r) for λ in tl, r in tr]
        @test size(interp(bounds, radii, table, tl, tr)) == (3, 2)
    end

    @testset "input order does not matter" begin
        tl, tr = [0.33, 0.7], [0.6, 1.7]
        ref = interp(bounds, radii, table, tl, tr)
        @test interp(reverse(bounds), radii, table[end:-1:1, :], tl, tr) ≈ ref
        @test interp(bounds, reverse(radii), table[:, end:-1:1], tl, tr) ≈ ref
        @test interp([collect(b) for b in bounds], radii, table, tl, tr) ≈ ref
    end

    @testset "extrapolation" begin
        tl, tr = [0.05, 10.0], [0.1, 20.0]  # Outside the table on both sides.
        exact = [powerlaw(λ, r) for λ in tl, r in tr]
        @test interp(bounds, radii, table, tl, tr) ≈ exact
        @test interp(bounds, radii, table, tl, tr; extrapolation=:linear) ≈ exact
        flat = interp(bounds, radii, table, tl, tr; extrapolation=:flat)
        @test flat ≈ [table[1, 1] table[1, end]; table[end, 1] table[end, end]]
        @test_throws BoundsError interp(bounds, radii, table, tl, tr; extrapolation=:error)
        @test_throws ArgumentError interp(
            bounds, radii, table, tl, tr; extrapolation=:cubic
        )
    end

    @testset "zeros and floor" begin
        z = copy(table)
        z[1, :] .= 0
        result = interp(bounds, radii, z, centers, radii)
        @test all(iszero, result[1, :])
        @test result[2:end, :] ≈ table[2:end, :]
        @test all(iszero, interp(bounds, radii, zero(table), [0.7], [1.0]))
        r32 = interp(bounds, Float32.(radii), Float32.(table), [0.7f0], [1.0f0])
        @test r32[1, 1] ≈ powerlaw(0.7, 1.0) rtol = 1e-5
    end

    @testset "invalid input" begin
        @test_throws ArgumentError interp(bounds[1:1], radii, table[1:1, :], [0.5], [1.0])
        @test_throws ArgumentError interp(bounds, radii[1:1], table[:, 1:1], [0.5], [1.0])
        @test_throws ArgumentError interp(bounds, radii, table[1:3, :], [0.5], [1.0])
        bad = [(0.0, 0.4); bounds[2:end]]
        @test_throws ArgumentError interp(bad, radii, table, [0.5], [1.0])
        bad = [(0.4, 0.2); bounds[2:end]]
        @test_throws ArgumentError interp(bad, radii, table, [0.5], [1.0])
        bad = [bounds[1:3]; [bounds[3]]]  # Repeated band.
        @test_throws ArgumentError interp(bad, radii, table, [0.5], [1.0])
        @test_throws ArgumentError interp(bounds, [0.5, 1.0, 1.0], table, [0.5], [1.0])
        @test_throws ArgumentError interp(bounds, [-0.5, 1.0, 2.0], table, [0.5], [1.0])
        @test_throws ArgumentError interp(bounds, radii, table, [0.0], [1.0])
        @test_throws ArgumentError interp(bounds, radii, table, [0.5], [-1.0])
        @test_throws ArgumentError interp(bounds, radii, table, [0.5], [1.0]; eps_val=0)
    end
end
