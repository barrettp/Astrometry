#### Astronomy / Galactic Coordinates

"""
    g2icrs(lon::Float64, lat::Float64)

Transformation from Galactic Coordinates to ICRS.

# Input

 - `lon`   -- galactic longitude (radians)
 - `lat`   -- galactic latitude (radians)

# Output

 - `ras`   -- ICRS right ascension (radians)
 - `dec`   -- ICRS declination (radians)

# Note

1) The IAU 1958 system of Galactic coordinates was defined with
   respect to the now obsolete reference system FK4 B1950.0.  When
   interpreting the system in a modern context, several factors have
   to be taken into account:

   . The inclusion in FK4 positions of the E-terms of aberration.

   . The distortion of the FK4 proper motion system by differential
     Galactic rotation.

   . The use of the B1950.0 equinox rather than the now-standard
     J2000.0.

   . The frame bias between ICRS and the J2000.0 mean place system.

   The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
   matrix that transforms directly between ICRS and Galactic
   coordinates with the above factors taken into account.  The matrix
   is derived from three angles, namely the ICRS coordinates of the
   Galactic pole and the longitude of the ascending node of the
   galactic equator on the ICRS equator.  They are given in degrees to
   five decimal places and for canonical purposes are regarded as
   exact.  In the Hipparcos Catalogue the matrix elements are given to
   10 decimal places (about 20 microarcsec).  In the present ERFA
   function the matrix elements have been recomputed from the
   canonical three angles and are given to 30 decimal places.

2) The inverse transformation is performed by the function eraIcrs2g.

# References

Perryman M.A.C. & ESA, 1997, ESA SP-1200, The Hipparcos and Tycho
catalogues.  Astrometric and photometric star catalogues derived from
the ESA Hipparcos Space Astrometry Mission.  ESA Publications
Division, Noordwijk, Netherlands.
"""
function g2icrs(lon::Float64, lat::Float64)
    #=
    L2,B2 system of galactic coordinates in the form presented in the
    Hipparcos Catalogue.  In degrees:
  
    P = 192.85948    right ascension of the Galactic north pole in ICRS
    Q =  27.12825    declination of the Galactic north pole in ICRS
    R =  32.93192    Galactic longitude of the ascending node of
                     the Galactic equator on the ICRS equator
  
    ICRS to galactic rotation matrix, obtained by computing Rz(-R)
    Rx(π/2-Q) Rz(π/2+P) to the full precision shown.
    =#
    ras, dec = c2s(r_gal_icrs'*s2c(lon, lat))
    NamedTuple{(:ras, :dec)}((anp(ras), anpm(dec)))
end

"""
    icrs2g(ras::Float64, dec::Float64)

Transformation from ICRS to Galactic Coordinates.

# Input

 - `ras`   -- ICRS right ascension (radians)
 - `dec`   -- ICRS declination (radians)

# Output

 - `lon`   -- galactic longitude (radians)
 - `lat`   -- galactic latitude (radians)

# Note

1) The IAU 1958 system of Galactic coordinates was defined with
   respect to the now obsolete reference system FK4 B1950.0.  When
   interpreting the system in a modern context, several factors have
   to be taken into account:

   . The inclusion in FK4 positions of the E-terms of aberration.

   . The distortion of the FK4 proper motion system by differential
     Galactic rotation.

   . The use of the B1950.0 equinox rather than the now-standard
     J2000.0.

   . The frame bias between ICRS and the J2000.0 mean place system.

   The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
   matrix that transforms directly between ICRS and Galactic
   coordinates with the above factors taken into account.  The matrix
   is derived from three angles, namely the ICRS coordinates of the
   Galactic pole and the longitude of the ascending node of the
   galactic equator on the ICRS equator.  They are given in degrees to
   five decimal places and for canonical purposes are regarded as
   exact.  In the Hipparcos Catalogue the matrix elements are given to
   10 decimal places (about 20 microarcsec).  In the present ERFA
   function the matrix elements have been recomputed from the
   canonical three angles and are given to 30 decimal places.

2) The inverse transformation is performed by the function eraG2icrs.

# References

Perryman M.A.C. & ESA, 1997, ESA SP-1200, The Hipparcos and Tycho
catalogues.  Astrometric and photometric star catalogues derived from
the ESA Hipparcos Space Astrometry Mission.  ESA Publications
Division, Noordwijk, Netherlands.
"""
function icrs2g(ras::Float64, dec::Float64)
    #=
    L2,B2 system of galactic coordinates in the form presented in the
    Hipparcos Catalogue.  In degrees:
  
    P = 192.85948    right ascension of the Galactic north pole in ICRS
    Q =  27.12825    declination of the Galactic north pole in ICRS
    R =  32.93192    Galactic longitude of the ascending node of
                     the Galactic equator on the ICRS equator
  
    ICRS to galactic rotation matrix, obtained by computing R_3(-R)
    R_1(π/2-Q) R_3(π/2+P) to the full precision shown:
    =#
    @inline lon, lat = c2s(r_gal_icrs*s2c(ras, dec))
    NamedTuple{(:lon, :lat)}((anp(lon), anpm(lat)))
end
