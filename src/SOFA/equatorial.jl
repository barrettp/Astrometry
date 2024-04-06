#### Astronomy / Horizontal-Equatorial

"""
    ae2hd(azimuth::Float64, altitude::Float64, latitude::Float64)

Horizon to equatorial coordinates. Transform azimuth and altitude to
hour angle and declination

# Input

 - `azimuth`  -- azimuth
 - `altitude` -- altitude (informally, elevation)
 - `latitude` -- site latitude

# Output

 - `equatorial` -- hour angle and declination

# Note

1) All the arguments are angles in radians.

2) The sign convention for azimuth is north zero, east +pi/2.

3) HA is returned in the range +/-π.  Declination is returned in the
   range +/-π/2.

4) The latitude ϕ is π/2 minus the angle between the Earth's rotation
   axis and the adopted zenith.  In many applications it will be
   sufficient to use the published geodetic latitude of the site.  In
   very precise (sub-arcsecond) applications, ϕ can be corrected for
   polar motion.

5) The azimuth az must be with respect to the rotational north pole,
   as opposed to the ITRS pole, and an azimuth with respect to north
   on a map of the Earth's surface will need to be adjusted for polar
   motion if sub-arcsecond accuracy is required.

6) Should the user wish to work with respect to the astronomical
   zenith rather than the geodetic zenith, ϕ will need to be adjusted
   for deflection of the vertical (often tens of arcseconds), and the
   zero point of ha will also be affected.

7) The transformation is the same as Ve = Ry(ϕ-π/2)*Rz(π)*Vh, where Ve
   and Vh are lefthanded unit vectors in the (ha,dec) and (az,el)
   systems respectively and Rz and Ry are rotations about first the
   z-axis and then the y-axis.  (n.b. Rz(π) simply reverses the signs
   of the x and y components.)  For efficiency, the algorithm is
   written out rather than calling other utility functions.  For
   applications that require even greater efficiency, additional
   savings are possible if constant terms such as functions of
   latitude are computed once and for all.

8) Again for efficiency, no range checking of arguments is carried out.
"""
function ae2hd(azimuth::Float64, altitude::Float64, latitude::Float64)
    sa, ca = sincos(azimuth)
    se, ce = sincos(altitude)
    sp, cp = sincos(latitude)

    x, y, z = -ca*ce*sp + se*cp, -sa*ce, ca*ce*cp + se*sp
    r = sqrt(x*x + y*y)
    (r != 0.0 ? atan(y, x) : 0.0, atan(z, r))
end

"""
    hd2ae(HA::Float64, Dec::Float64, ϕ::Float64)

Equatorial to horizon coordinates: transform hour angle and
declination to azimuth and altitude.

# Input

 - `HA`    -- hour angle (local)
 - `Dec`   -- declination
 - `ϕ`     -- site latitude

# Output

 - `azimuth`  -- azimuth
 - `altitude` -- altitude (informally, elevation)

# Note

1) All the arguments are angles in radians.

2) Azimuth is returned in the range 0-2pi; north is zero, and east is
   +pi/2.  Altitude is returned in the range +/- pi/2.

3) The latitude ϕ is pi/2 minus the angle between the Earth's rotation
   axis and the adopted zenith.  In many applications it will be
   sufficient to use the published geodetic latitude of the site.  In
   very precise (sub-arcsecond) applications, ϕ can be corrected for
   polar motion.

4) The returned azimuth az is with respect to the rotational north
   pole, as opposed to the ITRS pole, and for sub-arcsecond accuracy
   will need to be adjusted for polar motion if it is to be with
   respect to north on a map of the Earth's surface.

5) Should the user wish to work with respect to the astronomical
   zenith rather than the geodetic zenith, ϕ will need to be adjusted
   for deflection of the vertical (often tens of arcseconds), and the
   zero point of the hour angle ha will also be affected.

6) The transformation is the same as Vh = Rz(pi)*Ry(pi/2-ϕ)*Ve, where
   Vh and Ve are lefthanded unit vectors in the (az,el) and (ha,dec)
   systems respectively and Ry and Rz are rotations about first the
   y-axis and then the z-axis.  (n.b. Rz(pi) simply reverses the signs
   of the x and y components.)  For efficiency, the algorithm is
   written out rather than calling other utility functions.  For
   applications that require even greater efficiency, additional
   savings are possible if constant terms such as functions of
   latitude are computed once and for all.

7) Again for efficiency, no range checking of arguments is carried out.
"""
function hd2ae(HA::Float64, Dec::Float64, ϕ::Float64)
    sh, ch = sincos(HA)
    sd, cd = sincos(Dec)
    sp, cp = sincos(ϕ)

    x, y, z = -ch*cd*sp + sd*cp, -sh*cd, ch*cd*cp + sd*sp
    r = sqrt(x*x + y*y)
    (azi = mod2pi(r != 0.0 ? atan(y, x) : 0.0), alt = atan(z, r))
end

"""
    hd2pa(HA::Float64, Dec::Float64, latitude::Float64)

Parallactic angle for a given hour angle and declination.

# Input

 - `HA`    -- hour angle
 - `Dec`   -- declination
 - `latitude` -- site latitude

# Output

 - `angle` -- parallactic angle

# Note

1) All the arguments are angles in radians.

2) The parallactic angle at a point in the sky is the position angle
   of the vertical, i.e. the angle between the directions to the north
   celestial pole and to the zenith respectively.

3) The result is returned in the range -pi to +pi.

4) At the pole itself a zero result is returned.

5) The latitude ϕ is pi/2 minus the angle between the Earth's rotation
   axis and the adopted zenith.  In many applications it will be
   sufficient to use the published geodetic latitude of the site.  In
   very precise (sub-arcsecond) applications, ϕ can be corrected for
   polar motion.

6) Should the user wish to work with respect to the astronomical
   zenith rather than the geodetic zenith, ϕ will need to be adjusted
   for deflection of the vertical (often tens of arcseconds), and the
   zero point of the hour angle ha will also be affected.

# References

Smart, W.M., "Spherical Astronomy", Cambridge University Press, 6th
edition (Green, 1977), p49.
"""
function hd2pa(HA::Float64, Dec::Float64, latitude::Float64)
    sqsz = cos(latitude)*sin(HA)
    cqsz = sin(latitude)*cos(Dec) - cos(latitude)*sin(Dec)*cos(HA)
    sqsz != 0.0 || cqsz || 0.0 ? atan(sqsz, cqsz) : 0.0
end
