#### Astronomy / Star Catalogs

"""
    fk425(r1950::Float64, d1950::Float64, dr1950::Float64, dd1950::Float64,
          p1950::Float64, v1950::Float64)

Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

This function converts a star's catalog data from the old FK4
(Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system.

# Input (all B1950.0, FK4)

 - `r1950`  -- B1950.0 RA (rad)
 - `d1950`  -- B1950.0 Dec (rad)
 - `dr1950` -- B1950.0 RA proper moation (rad/trop-year)
 - `dd1950` -- B1950.0 Dec proper motions (rad/trop-year)
 - `p1950`  -- parallax (arcsec)
 - `v1950`  -- radial velocity (km/s, +ve = moving away)

# Output (all J2000.0, FK5)

 - `r2000`  -- J2000.0 RA (rad)
 - `d2000`  -- J2000.0 Dec (rad)
 - `dr2000` -- J2000.0 RA proper motion (rad/Jul-year)
 - `dd2000` -- J2000.0 Dec proper motions (rad/Jul-year)
 - `p2000`  -- parallax (arcsec)
 - `v2000`  -- radial velocity (km/s, +ve = moving away)

# Note

1) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
   and are per year rather than per century.

2) The conversion is somewhat complicated, for several reasons:

   . Change of standard epoch from B1950.0 to J2000.0.

   . An intermediate transition date of 1984 January 1.0 TT.

   . A change of precession model.

   . Change of time unit for proper motion (tropical to Julian).

   . FK4 positions include the E-terms of aberration, to simplify
     the hand computation of annual aberration.  FK5 positions
     assume a rigorous aberration computation based on the Earth's
     barycentric velocity.

   . The E-terms also affect proper motions, and in particular cause
     objects at large distances to exhibit fictitious proper
     motions.

   The algorithm is based on Smith et al. (1989) and Yallop et al.
   (1989), which presented a matrix method due to Standish (1982) as
   developed by Aoki et al. (1983), using Kinoshita's development of
   Andoyer's post-Newcomb precession.  The numerical constants from
   Seidelmann (1992) are used canonically.

3) Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
   Conversions for different epochs and equinoxes would require
   additional treatment for precession, proper motion and E-terms.

4) In the FK4 catalog the proper motions of stars within 10 degrees of
   the poles do not embody differential E-terms effects and should,
   strictly speaking, be handled in a different manner from stars
   outside these regions.  However, given the general lack of
   homogeneity of the star data available for routine astrometry, the
   difficulties of handling positions that may have been determined
   from astrometric fields spanning the polar and non- polar regions,
   the likelihood that the differential E-terms effect was not taken
   into account when allowing for proper motion in past astrometry,
   and the undesirability of a discontinuity in the algorithm, the
   decision has been made in this ERFA algorithm to include the
   effects of differential E-terms on the proper motions for all
   stars, whether polar or not.  At epoch J2000.0, and measuring "on
   the sky" rather than in terms of RA change, the errors resulting
   from this simplification are less than 1 milliarcsecond in position
   and 1 milliarcsecond per century in proper motion.

# References

Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0 FK4-based
positions of stars to epoch J2000.0 positions in accordance with the
new IAU resolutions".  Astron.Astrophys.  128, 263-267.

Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
Astronomical Almanac", ISBN 0-935702-68-7.

Smith, C.A. et al., 1989, "The transformation of astrometric catalog
systems to the equinox J2000.0".  Astron.J. 97, 265.

Standish, E.M., 1982, "Conversion of positions and proper motions from
B1950.0 to the IAU system at J2000.0".  Astron.Astrophys., 115, 1,
20-22.

Yallop, B.D. et al., 1989, "Transformation of mean star places from
FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".  Astron.J. 97,
274.
"""
function fk425(r1950::Float64, d1950::Float64, dr1950::Float64, dd1950::Float64,
	       p1950::Float64, v1950::Float64)
    #  Canonical constants (Seidelmann 1992)
    #  km/s to AU per tropical century.
    PMF, VF = 100.0*3600.0*rad2deg(1.0), 21.095
    #  constant pv-vector (cf. Seidelmann eq. 3.591-2, vectors A and Adot)
    a = [[-1.62557e-6, -0.31919e-6, -0.13843e-6],
         [+1.245e-3,   -1.580e-3,   -0.659e-3]]
    #  3x2 matrix of pv-vectors (cf. Seidelmann eq. 3.591-4, matrix M)
    em = [[[[+0.9999256782,     -0.0111820611,     -0.0048579477    ],
            [+0.00000242395018, -0.00000002710663, -0.00000001177656]],
           [[+0.0111820610,     +0.9999374784,     -0.0000271765    ],
            [+0.00000002710663, +0.00000242397878, -0.00000000006587]],
           [[+0.0048579479,     -0.0000271474,     +0.9999881997,   ],
            [+0.00000001177656, -0.00000000006582, +0.00000242410173]]],
          [[[-0.000551,         -0.238565,         +0.435739        ],
            [+0.99994704,       -0.01118251,       -0.00485767      ]],
           [[+0.238514,         -0.002667,         -0.008541        ],
            [+0.01118251,       +0.99995883,       -0.00002718      ]],
           [[-0.435623,         +0.012254,         +0.002117        ],
            [+0.00485767,       -0.00002714,       +1.00000956      ]]]]
    #  FK$ data (units radians and arcsec per tropical century).
    #  Express as a pv-vector
    r0 = s2pv(r1950, d1950, 1.0, PMF*dr1950, PMF*dd1950, VF*p1950*v1950)
    #  Allow for E-terms (cf. Seidelmann eq. 3.591-2).
    pv1 = pvmpv(r0, a)
    pv2 = [r0[1]'*a[1].*r0[1], r0[1]'*a[2].*r0[1]]
    pv1 = pvppv(pv1, pv2)
    #  Convert pv-vector to Fricke system (cf., Seidelmann eq. 3.591-3).
    pv2 = [em[1][1][1]'*pv1[1], em[1][1][2]'*pv2[2]]
end

function fk45z()
end

function fk524()
end

function fk52h()
end

function fk54z()
end

function fk5hip()
end

function fk5hz()
end

function h2fk5()
end

function hfk5z()
end

"""
    starpm(ras::Float64, dec::Float64, pmras::Float64, pmdec::Float64,
           plx::Float64, rvel::Float64, epoch1a::Float64, epoch1b::Float64,
           epoch2a::Float64, epoch2b::Float64)

Star proper motion:  update star catalog data for space motion.

# Input

 - `ra1`   -- right ascension (radians), before
 - `dec1`  -- declination (radians), before
 - `pmr1`  -- RA proper motion (radians/year), before
 - `pmd1`  -- Dec proper motion (radians/year), before
 - `px1`   -- parallax (arcseconds), before
 - `rv1`   -- radial velocity (km/s, +ve = receding), before
 - `ep1a`  -- "before" epoch, part A (Note 1)
 - `ep1b`  -- "before" epoch, part B (Note 1)
 - `ep2a`  -- "after" epoch, part A (Note 1)
 - `ep2b`  -- "after" epoch, part B (Note 1)

# Output

 - `ra2`   -- right ascension (radians), after
 - `dec2`  -- declination (radians), after
 - `pmr2`  -- RA proper motion (radians/year), after
 - `pmd2`  -- Dec proper motion (radians/year), after
 - `px2`   -- parallax (arcseconds), after
 - `rv2`   -- radial velocity (km/s, +ve = receding), after

   -1 = system error (should not occur)
    0 = no warnings or errors
    1 = distance overridden (Note 6)
    2 = excessive velocity (Note 7)
    4 = solution didn't converge (Note 8)
    else = binary logical OR of the above warnings

# Note

1) The starting and ending TDB dates ep1a+ep1b and ep2a+ep2b are
   Julian Dates, apportioned in any convenient way between the two
   parts (A and B).  For example, JD(TDB)=2450123.7 could be expressed
   in any of these ways, among others:

           epna          epnb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 method is best matched to the way the
   argument is handled internally and will deliver the optimum
   resolution.  The MJD method and the date & time methods are both
   good compromises between resolution and convenience.

2) In accordance with normal star-catalog conventions, the object's
   right ascension and declination are freed from the effects of
   secular aberration.  The frame, which is aligned to the catalog
   equator and equinox, is Lorentzian and centered on the SSB.

   The proper motions are the rate of change of the right ascension
   and declination at the catalog epoch and are in radians per TDB
   Julian year.

   The parallax and radial velocity are in the same frame.

3) Care is needed with units.  The star coordinates are in radians and
   the proper motions in radians per Julian year, but the parallax is
   in arcseconds.

4) The RA proper motion is in terms of coordinate angle, not true
   angle.  If the catalog uses arcseconds for both RA and Dec proper
   motions, the RA proper motion will need to be divided by cos(Dec)
   before use.

5) Straight-line motion at constant speed, in the inertial frame, is
   assumed.

6) An extremely small (or zero or negative) parallax is interpreted to
   mean that the object is on the "celestial sphere", the radius of
   which is an arbitrary (large) value (see the eraStarpv function for
   the value used).  When the distance is overridden in this way, the
   status, initially zero, has 1 added to it.

7) If the space velocity is a significant fraction of c (see the
   constant VMAX in the function eraStarpv), it is arbitrarily set to
   zero.  When this action occurs, 2 is added to the status.

8) The relativistic adjustment carried out in the eraStarpv function
   involves an iterative calculation.  If the process fails to
   converge within a set number of iterations, 4 is added to the
   status.
"""
function starpm(ras::Float64, dec::Float64, pmras::Float64, pmdec::Float64,
                plx::Float64, rvel::Float64, epoch1a::Float64, epoch1b::Float64,
                epoch2a::Float64, epoch2b::Float64)
    DC = 173.1446333113497
    #  Position-velocity vector
    pv1 = starpv(ras, dec, pmras, pmdec, plx, rvel)
    #  Light time when observed (days)
    tl1 = pm(pv1[1])/DC # (SECPERDAY*LIGHTSPEED/ASTRUNIT)
    #  Time interval
    Δt = (epoch2a - epoch1a) + (epoch2b - epoch1b)
    #  Move star along track from the "before" observed position to the
    #  "after" geometric position.
    pv = pvu(Δt + tl1, pv1)
    #  From this geometric position, deduce the observed light time (days)
    #  at the "after" epoch (with theoretically unnecesary error check).
    r2, rdv, v2 = pv[1]'*pv[1], pv[1]'*pv[2], pv[2]'*pv[2]
    c2mv2 = DC^2 - v2 # (SECPERDAY*LIGHTSPEED/ASTRUNIT)^2 - v2
    @assert c2mv2 > 0
    tl2 = (-rdv + sqrt(rdv^2 + c2mv2*r2))/c2mv2
    #  Move the position along track from the observed place at the
    #  "before" epoch to the observed place at the "after" epoch.
    pvstar(pvu(Δt + (tl1 - tl2), pv1))
end
