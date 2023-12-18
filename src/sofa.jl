"""
    Drift second parameters
"""
struct Driftsecond
    year::Integer
    month::Integer
    mjd::Float64
    offset::Float64
    rate::Float64
end

"""
    Leap second parameters
"""
struct Leapsecond
    year::Integer
    month::Integer
    second::Float64
end

"""
    Star independent astrometry parameters
"""
struct Astrom
    pmt::Float64          # proper motion time interval (SSB, Julian years)
    eb::Vector{Float64}  # SSB to observer (vector, AU)
    eh::Vector{Float64}  # Sun to observer (vector, unit)
    em::Float64          # distance from Sun to observer (AU)
    v::Vector{Float64}   # barycentric observer velocity (vector, c)
    bm1::Float64         # inverse Lorenz factor, i.e., sqrt(1-v^2) 
    bpn::Matrix{Float64} # bias-precesson-nutation matrix
    along::Float64       # longitude + s' + dERA(DUT) (radians)
    phi::Float64         # geodetic latitude (radians)
    xpl::Float64         # polar motion xp wrt local meridian (radians)
    ypl::Float64         # polar motion yp wrt local meridian (radians)
    sphi::Float64        # sine of geodetic latitude
    cphi::Float64        # cosine of geodetic latitude
    diurab::Float64      # magnitude of diurnal aberration vector
    eral::Float64        # `local` Earth rotation angle (radians)
    refa::Float64        # refraction constant A (radians)
    refb::Float64        # refraction constant B (radians)
end

function Astrom()
    Astrom(0., [0., 0., 0.], [0., 0., 0.], 0., [0., 0., 0.], 0.,
           zeros(Float64,3,3), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
end
function Astrom(pm, eb, eh, em, v, bm1, bpn)
    Astrom(pm, eb, eh, em, v, bm1, bpn,
           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
    Body parameters for light deflection
"""
struct Ldbody
    bm::Float64                   #  mass of the body (solar masses)
    dl::Float64                   #  deflection limiter (radians^2/2)
    pv::Vector{Vector{Float64}}   #  barycentric PV of the body (AU, AU/day)
end

const TINY = 1e-6

# Reference dates and drift rates (sec/day), pre leap seconds.
const DRIFTSECOND =
    [Driftsecond(1960,  1, 37300.0, 1.4178180, 0.0012960),
     Driftsecond(1961,  1, 37300.0, 1.4228180, 0.0012960),
     Driftsecond(1961,  8, 37300.0, 1.3728180, 0.0012960),
     Driftsecond(1962,  1, 37665.0, 1.8458580, 0.0011232),
     Driftsecond(1963, 11, 37665.0, 1.9458580, 0.0011232),
     Driftsecond(1964,  1, 38761.0, 3.2401300, 0.0012960),
     Driftsecond(1964,  4, 38761.0, 3.3401300, 0.0012960),
     Driftsecond(1964,  9, 38761.0, 3.4401300, 0.0012960),
     Driftsecond(1965,  1, 38761.0, 3.5401300, 0.0012960),
     Driftsecond(1965,  3, 38761.0, 3.6401300, 0.0012960),
     Driftsecond(1965,  7, 38761.0, 3.7401300, 0.0012960),
     Driftsecond(1965,  9, 38761.0, 3.8401300, 0.0012960),
     Driftsecond(1966,  1, 39126.0, 4.3131700, 0.0025920),
     Driftsecond(1968,  2, 39126.0, 4.2131700, 0.0025920)]

# Dates and Δ(AT)s.
const LEAPSECOND =
    [Leapsecond(1972,  1, 10.0),
     Leapsecond(1972,  7, 11.0),
     Leapsecond(1973,  1, 12.0),
     Leapsecond(1974,  1, 13.0),
     Leapsecond(1975,  1, 14.0),
     Leapsecond(1976,  1, 15.0),
     Leapsecond(1977,  1, 16.0),
     Leapsecond(1978,  1, 17.0),
     Leapsecond(1979,  1, 18.0),
     Leapsecond(1980,  1, 19.0),
     Leapsecond(1981,  7, 20.0),
     Leapsecond(1982,  7, 21.0),
     Leapsecond(1983,  7, 22.0),
     Leapsecond(1985,  7, 23.0),
     Leapsecond(1988,  1, 24.0),
     Leapsecond(1990,  1, 25.0),
     Leapsecond(1991,  1, 26.0),
     Leapsecond(1992,  7, 27.0),
     Leapsecond(1993,  7, 28.0),
     Leapsecond(1994,  7, 29.0),
     Leapsecond(1996,  1, 30.0),
     Leapsecond(1997,  7, 31.0),
     Leapsecond(1999,  1, 32.0),
     Leapsecond(2006,  1, 33.0),
     Leapsecond(2009,  1, 34.0),
     Leapsecond(2012,  7, 35.0),
     Leapsecond(2015,  7, 36.0),
     Leapsecond(2017,  1, 37.0)]

include("SOFA/calendars.jl")
include("SOFA/astrometry.jl")
include("SOFA/ephemerides.jl")
include("SOFA/coefficients.jl")
include("SOFA/precession.jl")
include("SOFA/rotations.jl")
include("SOFA/spacemotion.jl")
include("SOFA/starcatalogs.jl")
include("SOFA/ecliptic.jl")
include("SOFA/galactic.jl")
include("SOFA/geocentric.jl")
include("SOFA/timescales.jl")
include("SOFA/equatorial.jl")
include("SOFA/gnomonic.jl")
include("SOFA/vectorops.jl")
