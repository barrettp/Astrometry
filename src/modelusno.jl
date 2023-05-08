
# Numerical coefficients of ψa, ωa, and χa, along with ϵ0, the obliquity at
# J2000.0, are 4-angle formulation from Capitaine et al. (2003), eqs. (4),
# (37), and (39).
const ϵ0 = 84381.406 # arcseconds

# USNO precession model coefficients (from NOVAS)
const ψAcoeff = [0., 5038.481507, -1.0790069, -0.00114045,  0.000132851, -0.0000000951]
const ωAcoeff = [ϵ0,   -0.025754,  0.0512623, -0.00772503, -0.000000467,  0.0000003337]
const χAcoeff = [0.,   10.556403, -2.3814292, -0.00121197,  0.000170663, -0.0000000560]

"""
    precession(jd1, pos1, jd2)

Precesses equatorial rectangular coordinates from one epoch to another. One of the
two epochs must be J2000.0. The coordinates are referred to the mean dynamical
equator and equinox of the two respective epochs. From NOVAS 3.1 .

# Arguments
- `jd1::Real`: TDB Julian date of first epoch. See Note 1 below.
- `pos1::AbstractVector`: Position vector in geocentric rectangulare coordinates
  referred to the mean dynamical equator and equinox of the first epoch.
- `jd2::Real`: TDB Julian date of second epoch. See Note 1 below.

# Returns
- `pos2::AbstractVector`: Position vector in geocentric rectangulare coordinates
  referred to the mean dynamical equator and equinox of the second epoch.

# Notes
- 1. Either `tjd1` or `tjd2` must be 2451545.0 (J2000.0) TDB.
"""
function precession(tjd1, pos1, tjd2)

    @assert tjd1 == JD2000 || tjd2 == JD2000

    ϵ0a = deg2rad(1/3600)*ϵ0
    
    # Compute either the forward or reverse rotation matrix for precession
    if (tjd2 == JD2000)
        Δt  = (tjd1 - tjd2)/DAYPERCEN
        ψA, ωA, χA = deg2rad(1/3600).*[
            Polynomial(ψAcoeff)(Δt), Polynomial(ωAcoeff)(Δt), Polynomial(χAcoeff)(Δt)]
        pos2 = (Rz(ψA) * Rx(-ωA) * Rz(-χA) * Rx(ϵ0a))' * pos1
    else
        Δt  = (tjd2 - tjd1)/DAYPERCEN
        ψA, ωA, χA = deg2rad(1/3600).*[
            Polynomial(ψAcoeff)(Δt), Polynomial(ωAcoeff)(Δt), Polynomial(χAcoeff)(Δt)]
        pos2 = (Rz(ψA) * Rx(-ωA) * Rz(-χA) * Rx(ϵ0a)) * pos1
    end
    return pos2
end
