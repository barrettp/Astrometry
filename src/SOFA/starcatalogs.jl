#### Astronomy / Star Catalogs

"""
    fk425(ra::Float64, ded::Float64, δra::Float64, δdec::Float64, plx::Float64,
          rv::Float64)

Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

This function converts a star's catalog data from the old FK4
(Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system.

# Input (all B1950.0, FK4)

 - `ra`     -- B1950.0 RA (rad)
 - `dec`    -- B1950.0 Dec (rad)
 - `δra`    -- B1950.0 RA proper moation (rad/trop-year)
 - `δdec`   -- B1950.0 Dec proper motions (rad/trop-year)
 - `plx`    -- parallax (arcsec)
 - `rv`     -- radial velocity (km/s, +ve = moving away)

# Output (all J2000.0, FK5)

 - `ra`     -- J2000.0 RA (rad)
 - `dec`    -- J2000.0 Dec (rad)
 - `δra   ` -- J2000.0 RA proper motion (rad/Jul-year)
 - `δdec`   -- J2000.0 Dec proper motions (rad/Jul-year)
 - `plx`    -- parallax (arcsec)
 - `rv`     -- radial velocity (km/s, +ve = moving away)

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
function fk425(ra::Float64, dec::Float64, δra::Float64, δdec::Float64,
               plx::Float64, rv::Float64)
    ####  Canonical constants (Seidelmann 1992)
    #  km/s to AU/trop-century
    PMF, VF, TINY = 3.6e5*rad2deg(1.), 21.095, 1e-30
    #  FK4 data (units radians and arcsec per tropical century)
    #  Express as a pv-vector
    # r0 = vcat(s2pv(ra, dec, 1.0, PMF*δra, PMF*δdec, VF*plx*rv)...)
    r0_ = s2pv(ra, dec, 1.0, PMF*δra, PMF*δdec, VF*plx*rv)
    # println(r0)
    # println(r0_)
    #  Allow for E-terms (cf. Seidelmann eq. 3.591-2) and convert pv-vector to
    #  Fricke system (cf., Seidelmann eq. 3.591-3).
    # r01 = vcat(r0[1:3], r0[1:3])
    # println(r01)
    # r1  = r0 .- A_fk4_fk5 .+ (r01'*A_fk4_fk5)*r01
    # println(r1)
    r1_ = [r0_[1] .- A_fk4_fk5[1:3] .+ (r0_[1]'*A_fk4_fk5[1:3])*r0_[1],
           r0_[2] .- A_fk4_fk5[4:6] .+ (r0_[1]'*A_fk4_fk5[4:6])*r0_[1]]
    # println(r1_)
    # pv = M_fk4_fk5*(r0 .- A_fk4_fk5 .+ (r01'*A_fk4_fk5)*r01)
    pv_ = [[Mfk4fk5[i][j][1]'*r1_[1] + Mfk4fk5[i][j][2]'*r1_[2] for j=1:3] for i=1:2]
    # println(pv)
    # println(pv_)
    #  Revert to catalog form
    # cat = pv2s([pv[1:3], pv[4:6]])
    cat = pv2s(pv_)
    if plx > TINY
        plx, rv = plx/cat[3], cat[6]/(VF*plx)
    end
    NamedTuple{(:RA, :Dec, :δRA, :δDec, :plx, :rv)}(
        (mod2pi(cat[1]), cat[2], cat[4]/PMF, cat[5]/PMF, plx, rv))
end

"""
    fk45z(ra::Float64, dec::Float64, epoch::Float64)

Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
proper motion in the FK5 system.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

This function converts a star's catalog data from the old FK4
(Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system, in
such a way that the FK5 proper motion is zero.  Because such a star
has, in general, a non-zero proper motion in the FK4 system, the
function requires the epoch at which the position in the FK4 system
was determined.

# Input

 - `ra`    -- B1950.0 FK4 RA at epoch (rad)
 - `dec`   -- B1950.0 FK4 Dec at epoch (rad)
 - `epoch` -- Besselian epoch (e.g. 1979.3)

# Output

 - `ra`    -- J2000.0 FK5 RA (rad)
 - `dec`   -- J2000.0 FK5 Dec (rad)

# Note

1) The epoch bepoch is strictly speaking Besselian, but if a Julian
   epoch is supplied the result will be affected only to a negligible
   extent.

2) The method is from Appendix 2 of Aoki et al. (1983), but using the
   constants of Seidelmann (1992).  See the function iauFk425 for a
   general introduction to the FK4 to FK5 conversion.

3) Conversion from equinox B1950.0 FK4 to equinox J2000.0 FK5 only is
   provided for.  Conversions for different starting and/or ending
   epochs would require additional treatment for precession, proper
   motion and E-terms.

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
   decision has been made in this SOFA algorithm to include the
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
"""
function fk45z(ra::Float64, dec::Float64, epoch::Float64)
    PMF = 3.6e5*rad2deg(1.)
    #  Spherical coordinates to p-vector, adjust p-vector to give zero proper
    #  motion in FK5.
    r1, p = s2c(ra, dec),  A_fk4_fk5[1:3] .+ (epoch - 1950.0)/PMF.*A_fk4_fk5[4:6]
    #  Remove E-terms
    p  = r1 .- (p .- (r1'*p).*r1)
    #  Convert to Fricke system pv-vector (cf. Seidelmann 3.591-3)
    pv = [[Mfk4fk5[i][j][1]'*p for j=1:3] for i=1:2]
    #  Allow for fictitious proper motion
    pv = pvu((epj(epb2jd(epoch)...) - 2000.0)/PMF, pv)
    #  Revert to spherical coordinates
    ra, dec = c2s(pv[1])
    NamedTuple{(:RA, :Dec)}((mod2pi(ra), dec))
end

"""
    fk524(ra::Float64, dec::Float64, δra::Float64, δdec::Float64, plx::Float64,
          rv::Float64)

Convert J2000.0 FK5 star catalog data to B1950.0 FK4.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Input (all J2000.0, FK5)

 - `ra`    -- J2000.0 RA (rad)
 - `dec`   -- J2000.0 Dec (rad)
 - `δra`   -- J2000.0 RA proper motions (rad/Jul.yr)
 - `δdec`  -- J2000.0 Dec proper motions (rad/Jul.yr)
 - `plx`   -- parallax (arcsec)
 - `rv`    -- radial velocity (km/s, +ve = moving away)

# Output (all B1950.0, FK4)

 - `ra`    -- B1950.0 RA (rad)
 - `dec`   -- B1950.0 Dec (rad)
 - `δra`   -- B1950.0 RA proper motions (rad/trop.yr)
 - `δdec`  -- B1950.0 Dec proper motions (rad/trop.yr)
 - `plx`   -- parallax (arcsec)
 - `rv`    -- radial velocity (km/s, +ve = moving away)

# Note

1) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
   and are per year rather than per century.

2) The conversion is somewhat complicated, for several reasons:

   . Change of standard epoch from J2000.0 to B1950.0.

   . An intermediate transition date of 1984 January 1.0 TT.

   . A change of precession model.

   . Change of time unit for proper motion (Julian to tropical).

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
   decision has been made in this SOFA algorithm to include the
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
function fk524(ra::Float64, dec::Float64, δra::Float64, δdec::Float64, plx::Float64,
               rv::Float64)
    ####  Canonical constants (Seidelmann 1992)
    #  km/s to AU/trop-century
    PMF, VF, TINY = 3.6e5*rad2deg(1.), 21.095, 1e-30
    #  Express as a pv-vector
    r0 = s2pv(ra, dec, 1.0, PMF*δra, PMF*δdec, VF*plx*rv)
    #  Convert pv-vector to Bessel-Newcomb system (cf. Seidelmann 3.592-1)
    r1 = [[Mfk5fk4[i][j][1]'*r0[1] + Mfk5fk4[i][j][2]'*r0[2] for j=1:3] for i=1:2]
    #  Apply E-terms (equivalent to Seidelmann 3.592-3, one iteration)
    #  Direction
    p1 = r1[1] .+ norm(r1[1])*A_fk4_fk5[1:3] .- (r1[1]'*A_fk4_fk5[1:3])*r1[1]
    #  Direction
    pv1 = r1[1] .+ norm(p1)*A_fk4_fk5[1:3] .- (r1[1]'*A_fk4_fk5[1:3])*r1[1]
    pv2 = r1[2] .+ norm(p1)*A_fk4_fk5[4:6] .- (r1[1]'*A_fk4_fk5[4:6])*pv1
    #  Revert to catalog form
    cat = pv2s([pv1, pv2])
    if plx > TINY
        plx, rv = plx/cat[3], cat[6]/(plx*VF)
    end
    NamedTuple{(:RA, :Dec, :δRA, :δDec, :plx, :rv)}(
        (mod2pi(cat[1]), cat[2], cat[4]/PMF, cat[5]/PMF, plx, rv))
end
"""
    fk52h(ra::Float64, dec::Float64, δra::Float64, δdec::Float64, plx::Float64,
          rv::Float64)

Transform FK5 (J2000.0) star data into the Hipparcos system.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Input

 - `r5`    -- J2000.0 RA (radians)
 - `d5`    -- J2000.0 Dec (radians)
 - `dr5`   -- J2000.0 proper motion in RA (dRA/dt, rad/Jyear)
 - `dd5`   -- J2000.0 proper motion in Dec (dDec/dt, rad/Jyear)
 - `px5`   -- parallax (arcsec)
 - `rv5`   -- radial velocity (km/s, positive = receding)

Returned (all Hipparcos, epoch J2000.0):

 - `rh`    -- J2000.0 RA (Hipparcos; radians)
 - `dh`    -- J2000.0 Dec (Hipparcos; radians)
 - `drh`   -- J2000.0 proper motion in RA (dRA/dt, rad/Jyear)
 - `ddh`   -- J2000.0 proper motion in Dec (dDec/dt, rad/Jyear)
 - `pxh`   -- parallax (arcsec)
 - `rvh`   -- radial velocity (km/s, positive = receding)

# Note

1) This function transforms FK5 star positions and proper motions into
   the system of the Hipparcos catalog.

2) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
   and are per year rather than per century.

3) The FK5 to Hipparcos transformation is modeled as a pure rotation
   and spin; zonal errors in the FK5 catalog are not taken into
   account.

4) See also iauH2fk5, iauFk5hz, iauHfk5z.

Reference:

F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
"""
function fk52h(ra::Float64, dec::Float64, δra::Float64, δdec::Float64,
               plx::Float64, rv::Float64)
    #  FK5 to Hipparcos orientation matrix and rotation vector
    ϵ, ω = fk5hip()
    #  FK5 barycentric normalized pv-vector
    pv = starpv(ra, dec, δra, δdec, plx, rv)
    #  Orient the FK5 position into the Hipparcos system
    #  Apply spin to the position giving an extra space motion component, add it
    #  to the FK5 space motion, orient the FK5 space motion into the Hipparcos
    #  system, and convert Hipparcos pv-vector to spherical coordinates
    pvstar([ϵ*pv[1], ϵ*(vec2mat(pv[1])*ω/DAYPERYEAR .+ pv[2])])
end

"""
    fk54z(ra::Float64, dec::Float64, epoch::Float64)

Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
proper motion in FK5 and parallax.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Input

 - `ra`    -- J2000.0 FK5 RA (rad)
 - `dec`   -- J2000.0 FK5 Dec (rad)
 - `epoch` -- Besselian epoch (e.g. 1950.0)

# Output

 - `ra`    -- B1950.0 FK4 RA (rad) at epoch BEPOCH
 - `dec`   -- B1950.0 FK4 Dec (rad) at epoch BEPOCH
 - `δra`   -- B1950.0 FK4 proper motions in RA (rad/trop.yr)
 - `δdec`  -- B1950.0 FK4 proper motions in Dec (rad/trop.yr)

# Note

1) In contrast to the iauFk524 function, here the FK5 proper motions,
   the parallax and the radial velocity are presumed zero.

2) This function converts a star position from the IAU 1976 FK5
   (Fricke) system to the former FK4 (Bessel-Newcomb) system, for
   cases such as distant radio sources where it is presumed there is
   zero parallax and no proper motion.  Because of the E-terms of
   aberration, such objects have (in general) non-zero proper motion
   in FK4, and the present function returns those fictitious proper
   motions.

3) Conversion from J2000.0 FK5 to B1950.0 FK4 only is provided for.
   Conversions involving other equinoxes would require additional
   treatment for precession.

4) The position returned by this function is in the B1950.0 FK4
   reference system but at Besselian epoch bepoch.  For comparison
   with catalogs the bepoch argument will frequently be 1950.0. (In
   this context the distinction between Besselian and Julian epoch is
   insignificant.)

5) The RA component of the returned (fictitious) proper motion is
   dRA/dt rather than cos(Dec)*dRA/dt.
"""
function fk54z(ra::Float64, dec::Float64, epoch::Float64)
    #  FK5 equinox J2000.0 to FK4 equinox B1950.0
    cat = fk524(ra, dec, 0., 0., 0., 0.)
    #  Spherical to Cartesian
    p1 = s2c(cat[1], cat[2])
    #  Fictitious proper motion (radians/year) and apply the motion
    p1 .+= (epoch - 1950.0)*(
        cat[3]*[-p1[2], p1[1], 0.] .+
        cat[4]*[-cos(cat[1])*sin(cat[2]), -sin(cat[1])*sin(cat[2]), cos(cat[2])])
    #  Cartesian to spherical
    ra, dec = c2s(p1)
    NamedTuple{(:RA, :Dec, :δRA, :δDec)}((mod2pi(ra), dec, cat[3], cat[4]))
end

"""
    fk5hip()

FK5 to Hipparcos rotation and spin.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Output

 - `rmat`  -- FK5 rotation matrix wrt Hipparcos (Note 2)
 - `rot`   -- FK5 spin vector wrt Hipparcos (Note 3)

# Note

1) This function models the FK5 to Hipparcos transformation as a pure
   rotation and spin; zonal errors in the FK5 catalog are not taken
   into account.

2) The r-matrix r5h operates in the sense:

         P_Hipparcos = r5h x P_FK5

   where P_FK5 is a p-vector in the FK5 frame, and P_Hipparcos is the
   equivalent Hipparcos p-vector.

3) The r-vector s5h represents the time derivative of the FK5 to
   Hipparcos rotation.  The units are radians per year (Julian, TDB).

# References

F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
"""
function fk5hip()
    #  FK5 wrt Hipparcos orientation and spin (radians, radians/year)
    NamedTuple{(:pos, :rot)}(
        (rv2m(deg2rad(1/3600.0).*ϵ_hipparcos), deg2rad(1/3600.0).*ω_hipparcos))
end

"""
    fk5hz(ra::Float64, dec::Float64, day1::Float64, day2::Float64)

Transform an FK5 (J2000.0) star position into the system of the
Hipparcos catalog, assuming zero Hipparcos proper motion.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Input

 - `ra`    -- FK5 RA (radians), equinox J2000.0, at date
 - `dec`   -- FK5 Dec (radians), equinox J2000.0, at date
 - `day1`  -- TDB date (Notes 1,2)
 - `day2`  -- TDB date (Notes 1,2)

# Output

 - `ra`    -- Hipparcos RA (radians)
 - `dec`   -- Hipparcos Dec (radians)

# Note

1) This function converts a star position from the FK5 system to the
   Hipparcos system, in such a way that the Hipparcos proper motion is
   zero.  Because such a star has, in general, a non-zero proper
   motion in the FK5 system, the function requires the date at which
   the position in the FK5 system was determined.

2) The TT date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

          date1          date2

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

3) The FK5 to Hipparcos transformation is modeled as a pure rotation
   and spin; zonal errors in the FK5 catalog are not taken into
   account.

4) The position returned by this function is in the Hipparcos
   reference system but at date date1+date2.

5) See also iauFk52h, iauH2fk5, iauHfk5z.

# References

F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
"""
function fk5hz(ra::Float64, dec::Float64, day1::Float64, day2::Float64)
    #  Interval from given date to fundamental epoch J2000.0 (Julian years)
    Δt = - ((day1 - JD2000) + day2)/DAYPERYEAR
    #  FK5 to Hipparcos orientation matrix and rotation vector
    ϵ, ω = fk5hip()
    #  FK5 barycentric position vector, accumulate Hipparcos wrt FK5 spin over
    #  the interval, de-rotate the vector's FK5 axes back to the date, rotate
    #  the vector into the Hipparcos system, and convert to spherical
    ra, dec = c2s(ϵ*rv2m(Δt.*ω)'*s2c(ra, dec))
    NamedTuple{(:RA, :Dec)}((mod2pi(ra), dec))
end

"""
    h2fk5(ra::Float64, dec::Float64, δra::Float64, δdec::Float64, plx::Float64,
          rv::Float64)

Transform Hipparcos star data into the FK5 (J2000.0) system.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Input

 - `ra`    -- J2000.0 RA (Hipparcos; radians)
 - `dec`   -- J2000.0 Dec (Hipparcos; radians)
 - `δra`   -- J2000.0 proper motion in RA (dRA/dt, rad/Jyear)
 - `δdec`  -- J2000.0 proper motion in Dec (dDec/dt, rad/Jyear)
 - `plz`   -- parallax (arcsec)
 - `rv`    -- radial velocity (km/s, positive = receding)

# Output

 - `ra`    -- J2000.0 FK5 RA (radians)
 - `dec`   -- J2000.0 FK5 Dec (radians)
 - `δra`   -- J2000.0 FK5 proper motion in RA (dRA/dt, rad/Jyear)
 - `δdec`  -- J2000.0 FK5 proper motion in Dec (dDec/dt, rad/Jyear)
 - `plx`   -- parallax (arcsec)
 - `rv`    -- radial velocity (km/s, positive = receding)

# Note

1) This function transforms Hipparcos star positions and proper
   motions into FK5 J2000.0.

2) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
   and are per year rather than per century.

3) The FK5 to Hipparcos transformation is modeled as a pure rotation
   and spin; zonal errors in the FK5 catalog are not taken into
   account.

4) See also iauFk52h, iauFk5hz, iauHfk5z.

# References

F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
"""
function h2fk5(ra::Float64, dec::Float64, δra::Float64, δdec::Float64,
               plx::Float64, rv::Float64)
    #  Hipparcos barycentric normalized pv-vector
    pv = starpv(ra, dec, δra, δdec, plx, rv)
    #  FK5 to Hipparcos orientation matrix and spin vector
    ϵ, ω = fk5hip()
    #  Scale spin to units per day, orient the spin to the Hipparcos system,
    #  apply spin to the position giving an extra space motion component,
    #  subtract this component from the Hipparcos space motion, de-orient
    #  the Hipparcos space motion into the FK5 system, and convert the FK5
    # pv-vector to spherical coordinates.
    pvstar([ϵ'*pv[1], ϵ'*(pv[2] .- vec2mat(pv[1])*ϵ*ω./DAYPERYEAR)])
end

"""
    hfk5z(ra::Float64, dec::Float64, day1::Float64, day2::Float64)

Transform a Hipparcos star position into FK5 J2000.0, assuming zero
Hipparcos proper motion.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Input

 - `ra`    -- Hipparcos RA (radians)
 - `dec`   -- Hipparcos Dec (radians)
 - `day1`  -- TDB date (Note 1)
 - `day2`  -- TDB date (Note 1)

Returned (all FK5, equinox J2000.0, date date1+date2):

 - `ra`    -- RA (radians)
 - `dec`   -- Dec (radians)
 - `δra`   -- RA proper motion (rad/year, Note 4)
 - `δdec`  -- Dec proper motion (rad/year, Note 4)

# Note

1) The TT date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

          date1          date2

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

2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3) The FK5 to Hipparcos transformation is modeled as a pure rotation
   and spin; zonal errors in the FK5 catalog are not taken into
   account.

4) It was the intention that Hipparcos should be a close approximation
   to an inertial frame, so that distant objects have zero proper
   motion; such objects have (in general) non-zero proper motion in
   FK5, and this function returns those fictitious proper motions.

5) The position returned by this function is in the FK5 J2000.0
   reference system but at date date1+date2.

6) See also iauFk52h, iauH2fk5, iauFk5hz.

# References

F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
"""
function hfk5z(ra::Float64, dec::Float64, day1::Float64, day2::Float64)
    #  Time interval from fundamental epoch J2000.0 to given date (Julian year)
    Δt = ((day1 - JD2000) + day2)/DAYPERYEAR
    #  FK5 to Hipparcos orientiation matrix and rotation vector
    ϵ, ω = fk5hip()
    #  FK5 pv-vector to spherical
    p  = s2c(ra, dec)
    #  Accumulate Hipparcos wrt to FK5 spin over interval, de-orient & de-spin
    #  the Hipparcos position into FK5 J2000.0, rotate the spin to the
    #  Hipparcos system, apply spin to the position giving the space motion,
    #  and de-orient and de-spin the Hipparcos space motion into the Fk5 J2000.0
    cat = pv2s([(ϵ*rv2m(Δt*ω))'*p, (ϵ*rv2m(Δt*ω))'*(vec2mat(ϵ*ω)*p)])
    NamedTuple{(:RA, :Dec, :δRA, :δDec)}((mod2pi(cat[1]), cat[2], cat[4], cat[5]))
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
