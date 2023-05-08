module Astrometry

export JD2000, MJD0, calendar_mjd
export iau_1980_nutation, iau_1980_obliquity, iau_1976_precession
export iau_2000_equinox_complement, iau_2000a_nutation, iau_2000b_nutation
export iau_2000_position, iau_2000_precession
export iau_2006a_nutation, iau_2006_cio_locator, iau_2006_cip_xy, iau_2006_obliquity
export iau_2006_precession, iau_2006_tdb_tt
export equinox, precession_nutation, proper_motion, radec2rad

using Polynomials

include("constants.jl")
include("util.jl")
include("model1980.jl")
include("model2000.jl")
include("model2006.jl")
include("earth.jl")

end

#   Issues
#   1. Timing accuracy, i.e. significant digits.
