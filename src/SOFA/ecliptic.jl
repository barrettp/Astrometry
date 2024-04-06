#### Astronomy / Ecliptic Coordinates

"""
    eceq06(day1::Float64, day2::Float64, lon::Float64, lat::Float64)

Transformation from ecliptic coordinates (mean equinox and ecliptic of
date) to ICRS RA,Dec, using the IAU 2006 precession model.

# Input

 - `day1`  -- TT as a 2-part Julian date
 - `day2`  -- ... Julian date (Note 1)
 - `lon`   -- ecliptic longitude (radians)
 - `lat`   -- ecliptic latitude (radians)

# Output

 - `ras`   -- ICRS right ascension (radians)
 - `dec`   -- ICRS declination (radians)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

          day1          day2

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

2) No assumptions are made about whether the coordinates represent
   starlight and embody astrometric effects such as parallax or
   aberration.

3) The transformation is approximately that from ecliptic longitude
   and latitude (mean equinox and ecliptic of date) to mean J2000.0
   right ascension and declination, with only frame bias (always less
   than 25 mas) to disturb this classical picture.
"""
function eceq06(day1::Float64, day2::Float64, lon::Float64, lat::Float64)
    @inline ras, dec = c2s(ecm06(day1, day2)'*s2c(lon, lat))
    NamedTuple{(:ras, :dec)}((anp(ras), anpm(dec)))
end

"""
    ecm06(day1::Float64, day2::Float64)

ICRS equatorial to ecliptic rotation matrix, IAU 2006.

# Input

 - `day1`  -- TT as a 2-part Julian date (Note 1)
 - `day2`  -- ...Julian date (Note 1)

# Output

 - `r`     -- ICRS to ecliptic rotation matrix

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

          day1          day2

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

2) The matrix is in the sense

      E_ep = rm x P_ICRS,

   where P_ICRS is a vector with respect to ICRS right ascension and
   declination axes and E_ep is the same vector with respect to the
   (inertial) ecliptic and equinox of date.

3) P_ICRS is a free vector, merely a direction, typically of unit
   magnitude, and not bound to any particular spatial origin, such as
   the Earth, Sun or SSB.  No assumptions are made about whether it
   represents starlight and embodies astrometric effects such as
   parallax or aberration.  The transformation is approximately that
   between mean J2000.0 right ascension and declination and ecliptic
   longitude and latitude, with only frame bias (always less than 25
   mas) to disturb this classical picture.
"""
function ecm06(day1::Float64, day2::Float64)
    @inline Rx(obl06(day1, day2))*pmat06(day1, day2)
end

"""
    eqec06(day1::Float64, day2::Float64, ras::Float64, dec::Float64)
    
Transformation from ICRS equatorial coordinates to ecliptic
coordinates (mean equinox and ecliptic of date) using IAU 2006
precession model.

# Input

 - `day1`  -- TT as a 2-part Julian date (Note 1)
 - `day2`  -- ... Julian date (Note 1)
 - `ras`   -- ICRS right ascension (radians)
 - `dec`   -- ICRS declination (radians)

# Output

- `lon`    -- ecliptic longitude (radians)
- `lat`    -- ecliptic latitude (radians)

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

2) No assumptions are made about whether the coordinates represent
   starlight and embody astrometric effects such as parallax or
   aberration.

3) The transformation is approximately that from mean J2000.0 right
   ascension and declination to ecliptic longitude and latitude (mean
   equinox and ecliptic of date), with only frame bias (always less
   than 25 mas) to disturb this classical picture.
"""
function eqec06(day1::Float64, day2::Float64, ras::Float64, dec::Float64)
    @inline lon, lat = c2s(ecm06(day1, day2)*s2c(ras, dec))
    (lon = anp(lon), lat = anpm(lat))
end

"""
    lteceq(epoch::Float64, lon::Float64, lat::Float64)

Transformation from ecliptic coordinates (mean equinox and ecliptic of
date) to ICRS RA,Dec, using a long-term precession model.

# Input

 - `epoch` -- Julian epoch (TT)
 - `lon`   -- ecliptic longitude (radians)
 - `lat`   -- ecliptic latitude (radians)

# Output

 - `ras`   -- ICRS right ascension (radians)
 - `dec`   -- ICRS declination (radians)

# Note

1) No assumptions are made about whether the coordinates represent
   starlight and embody astrometric effects such as parallax or
   aberration.

2) The transformation is approximately that from ecliptic longitude
   and latitude (mean equinox and ecliptic of date) to mean J2000.0
   right ascension and declination, with only frame bias (always less
   than 25 mas) to disturb this classical picture.

3) The Vondrak et al. (2011, 2012) 400 millennia precession model
   agrees with the IAU 2006 precession at J2000.0 and stays within 100
   microarcseconds during the 20th and 21st centuries.  It is accurate
   to a few arcseconds throughout the historical period, worsening to
   a few tenths of a degree at the end of the +/- 200,000 year time
   span.

# References

Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
expressions, valid for long time intervals, Astron.Astrophys. 534, A22

Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
expressions, valid for long time intervals (Corrigendum),
Astron.Astrophys. 541, C1
"""
function lteceq(epoch::Float64, lon::Float64, lat::Float64)
    @inline ras, dec = c2s(ltecm(epoch)'*s2c(lon, lat))
    (RA = anp(ras), Dec = anpm(dec))
end

"""
    ltecm(epoch::Float64)

ICRS equatorial to ecliptic rotation matrix, long-term.

# Input

 - `epoch` -- Julian epoch (TT)

# Output

 - `r`     -- ICRS to ecliptic rotation matrix

# Note

1) The matrix is in the sense

      E_ep = rm x P_ICRS,

   where P_ICRS is a vector with respect to ICRS right ascension and
   declination axes and E_ep is the same vector with respect to the
   (inertial) ecliptic and equinox of epoch epj.

2) P_ICRS is a free vector, merely a direction, typically of unit
   magnitude, and not bound to any particular spatial origin, such as
   the Earth, Sun or SSB.  No assumptions are made about whether it
   represents starlight and embodies astrometric effects such as
   parallax or aberration.  The transformation is approximately that
   between mean J2000.0 right ascension and declination and ecliptic
   longitude and latitude, with only frame bias (always less than 25
   mas) to disturb this classical picture.

3) The Vondrak et al. (2011, 2012) 400 millennia precession model
   agrees with the IAU 2006 precession at J2000.0 and stays within 100
   microarcseconds during the 20th and 21st centuries.  It is accurate
   to a few arcseconds throughout the historical period, worsening to
   a few tenths of a degree at the end of the +/- 200,000 year time
   span.

# References

Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
expressions, valid for long time intervals, Astron.Astrophys. 534, A22

Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
expressions, valid for long time intervals (Corrigendum),
Astron.Astrophys. 541, C1
"""
function ltecm(epoch::Float64)
    #  Bias vector
    bias = [η0_2010, -ϵ0_2010, -dα0_2010]

    #  Equatorial and ecliptic poles
    @inline equ, ecl = ltpequ(epoch), ltpecl(epoch)
    
    #  Create matrix
    eqx = vec2mat(equ)*ecl/norm(vec2mat(equ)*ecl)
    vcat(eqx', (vec2mat(ecl)*eqx)', ecl')*(I + deg2rad(1/3600)*vec2mat(bias))
end

"""
    lteqec(epoch::Float64, ras::Float64, dec::Float64)

Transformation from ICRS equatorial coordinates to ecliptic
coordinates (mean equinox and ecliptic of date) using a long-term
precession model.

# Input

 - `epoch` -- Julian epoch (TT)
 - `ras`   -- ICRS right ascension (radians)
 - `dec`   -- ICRS declination (radians)

# Output

 - `lon`   -- ecliptic longitude (radians)
 - `lat`   -- ecliptic latitude (radians)

# Note

1) No assumptions are made about whether the coordinates represent
   starlight and embody astrometric effects such as parallax or
   aberration.

2) The transformation is approximately that from mean J2000.0 right
   ascension and declination to ecliptic longitude and latitude (mean
   equinox and ecliptic of date), with only frame bias (always less
   than 25 mas) to disturb this classical picture.

3) The Vondrak et al. (2011, 2012) 400 millennia precession model
   agrees with the IAU 2006 precession at J2000.0 and stays within 100
   microarcseconds during the 20th and 21st centuries.  It is accurate
   to a few arcseconds throughout the historical period, worsening to
   a few tenths of a degree at the end of the +/- 200,000 year time
   span.

# References

Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
expressions, valid for long time intervals, Astron.Astrophys. 534, A22

Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
expressions, valid for long time intervals (Corrigendum),
Astron.Astrophys. 541, C1
"""
function lteqec(epoch::Float64, ras::Float64, dec::Float64)
    @inline lon, lat = c2s(ltecm(epoch)*s2c(ras, dec))
    (lon = anp(lon), lat = anpm(lat))
end
