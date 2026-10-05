module PlanetAerosolUtils

export calc_size_integrated_effective_radius, interpolate_extinction_coeff_log

using Interpolations

"""
    calc_size_integrated_effective_radius(mass, density, radius)

Calculates the size-integrated effective radius (r_eff) as the cross-sectional
weighted average across multiple size bins.

This calculation is based on the ratio of the total volume of the aerosol
population to its total projected cross-sectional area.

# Arguments
- `mass::AbstractVector`: Mass of aerosol in each size bin [kg or g].
- `density::AbstractVector`: Particle density for each size bin [kg/m^3 or g/cm^3].
- `radius::AbstractVector`: Representative radius (e.g., geometric mean) of
  particles in each bin [m or microns].

# Returns
- A scalar value representing the effective radius of the total population in the
  same units as the input `radius`.

# Method
1. Calculate volume (V) for each bin: V = mass / density.
2. The effective radius is derived from the ratio of total volume to total
   area (V/A). For spherical particles, this simplifies to:
   r_eff = sum(V) / sum(V / r)

# Example
```julia
m = [0.1, 0.5, 0.2]
rho = [2000.0, 2000.0, 2000.0]
r = [0.1, 0.5, 2.0]
reff = calc_size_integrated_effective_radius(m, rho, r)
"""
function calc_size_integrated_effective_radius(mass::AbstractVector, density::AbstractVector,
        radius::AbstractVector)

    volume = mass ./ density
    numerator = sum(volume)
    denominator = sum(volume ./ radius)
    sizeIntegratedEffectiveRadius = numerator / denominator

    return sizeIntegratedEffectiveRadius
end # calc_size_integrated_effective_radius

"""
    interpolate_extinction_coeff_log(lambda_bounds, orig_radii, ext_matrix, target_lambda, target_radii; eps_val = 1e-12)

Interpolates mass extinction coefficients across a 2D grid of wavelengths and dust particle sizes using
log-log mapping. This method is specifically designed for aerosol optical properties, which often
vary exponentially across several orders of magnitude.

# Arguments
- `lambda_bounds::Vector{Tuple{Real, Real}}`: A vector of length `k` containing pairs of
  (lower_bound, upper_bound) for the original spectral bands. The function calculates the
  geometric mean (midpoint of the log values) for each band to determine its spectral center.
- `orig_radii::Vector{Real}`: A vector of length `l` containing the distinct dust effective
  radii (usually in micrometers) corresponding to the original lookup table.
- `ext_matrix::Matrix{Real}`: A `k × l` matrix where `ext_matrix[i, j]` is the mass extinction
  coefficient for the `i`-th wavelength band and `j`-th effective radius.
- `target_lambda::AbstractVector{Real}`: A vector of length `m` containing the target wavelengths
  to which the coefficients should be interpolated.
- `target_radii::AbstractVector{Real}`: A vector of length `n` containing the target dust effective
  radii to which the coefficients should be interpolated.
- `eps_val::Real`: A small positive floor value used to replace any zeros or negative values in
  the input matrix before log-transformation. Defaults to `1e-12`.

# Returns
- `new_ext_matrix::Matrix{Float64}`: An `m × n` matrix containing the interpolated mass
  extinction coefficients. The output is transformed back into linear space via the
  exponential function.
"""
function interpolate_extinction_coeff_log(
    lambda_bounds,
    orig_radii,
    ext_matrix,
    target_lambda,
    target_radii;
    eps_val = 1e-12
)
    # Coordinate Preparation
    orig_lambda_mid = [sqrt(b[1] * b[2]) for b in lambda_bounds]

    # Transform coordinates to Log-Space
    log_lambda_orig = log.(orig_lambda_mid)
    log_r_orig      = log.(orig_radii)

    # Transform data to Log-Space with safety floor
    log_ext_orig    = log.(max.(ext_matrix, eps_val))

    # Construct Interpolation Object with explicit Package prefix
    # Using Interpolations. prefix ensures it works inside your module
    itp = Interpolations.extrapolate(
        Interpolations.interpolate((log_lambda_orig, log_r_orig), log_ext_orig,
            Interpolations.Gridded(Interpolations.Linear())),
        Interpolations.Line()
    )

    # Evaluate Interpolation on the target grid
    new_ext_matrix = [
        exp(itp(log(max(tl, eps_val)), log(max(tr, eps_val))))
        for tl in target_lambda, tr in target_radii
    ]

    # Post-Processing
    new_ext_matrix[new_ext_matrix .<= eps_val * 1.01] .= 0.0

    return new_ext_matrix
end # function interpolate_extinction_coeff_log

end # module PlanetAerosolUtils
