#   IAU 2000 Model Common Constants

#   Initial obliquity of the ecliptic in arcseconds
const ϵ0_2000 = 84381.448

#   Precession model constants in mas (Lieske et al. 1977)
const ψ_1977 = [    0.0, 5038.7784, -1.07259, -0.001147]
const ω_1977 = [ϵ0_2000,    0.0,     0.05127, -0.007726]
const χ_1977 = [    0.0,   10.5526, -2.38064, -0.001125]

#   Frame bias model constants in mas 
const icrs_ra_2000 = -0.0146
const ψ_bias_2000  = -0.041775
const ϵ_bias_2000  = -0.0068192

#   Precession model correction constants in mas
const ψ_corr_2000 = -0.29965
const ϵ_corr_2000 = -0.02524

#   Matrix element for orienting the analytical model to DE405
const iau_2000_bcrs = [
    1.0             0.000000211284 -0.000000091603;
    -0.000000230286 0.917482137087 -0.397776982902;
    0.0             0.397776982902  0.917482137087]
