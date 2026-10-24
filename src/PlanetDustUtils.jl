module PlanetDustUtils

export calc_fractional_dust_masses_BFT

using SpecialFunctions

function calc_fractional_dust_masses_BFT(bounds)
    ## calculates fractional dust masses within an array of bounds
    ## based on the mass size distribution of brittle fragmentation theory
    ## (BFT, Kok, PNAS 2011, Eq. 6)

    # some parameters of BFT
    dAridSoils = 3.4
    sigmaAridSoils = 3.0
    Cv = 12.63 # (μm)
    lambda = 12 # (μm)

    nBins = length(bounds) - 1
    massFractions = Array{Float64}(undef, nBins)
    for i in 1:nBins
        bin_mean = sqrt(bounds[i] * bounds[i + 1])
        erf_in = log(bin_mean / dAridSoils) / (sqrt(2.0) * log(sigmaAridSoils))

        massFractions[i] =
            1.0 / Cv *
            (1 + erf(erf_in)) *
            exp(-(bin_mean / lambda)^3) *
            (bounds[i + 1] - bounds[i])
    end

    massFractions = massFractions / sum(massFractions)

    return massFractions
end # calc_fractional_dust_masses_BFT

end # module PlanetDustUtils
