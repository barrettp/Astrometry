####    Astronomy / Calendars    ####

"""
    cal2jd(year::Integer, month::Integer, day::Integer)

# Input

 - `year`  -- year in Gregorian calendar (Note 1)
 - `month` -- month in Gregorian calendar (Note 1)
 - `day`   -- day in Gregorian calendar (Note 1)

# Output

 - `mjd0`  -- MJD zero-point: always 2400000.5
 - `mjd`   -- Modified Julian Date for 0 hrs

# Note

1) The algorithm used is valid from -4800 March 1, but this
   implementation rejects dates before -4799 January 1.

2) The Julian Date is returned in two pieces, in the usual ERFA
   manner, which is designed to preserve time resolution.  The Julian
   Date is available as a single number by adding MJD0 and MJD.

3) In early eras the conversion is from the "Proleptic Gregorian
   Calendar"; no account is taken of the date(s) of adoption of the
   Gregorian Calendar, nor is the AD/BC numbering convention observed.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 12.92
(p604).
"""
function cal2jd(year::Integer, month::Integer, day::Integer)
    (mjd0 = MJDAY0, mjd = calendar2MJD(year, month, day))
end

"""
    epb(day1::Real, day2::Real)

Julian Date to Besselian Epoch.

# Input

 - `day1`   -- Julian Date (see note)
 - `day2`   -- Julian Date (see note)

# Output

 - `bessel` -- Besselian Epoch

# Note

The Julian Date is supplied in two pieces, in the usual ERFA manner,
which is designed to preserve time resolution.  The Julian Date is
available as a single number by adding day1 and day2.  The maximum
resolution is achieved if day1 is 2451545.0 (J2000.0).

# References

Lieske, J.H., 1979. Astron.Astrophys., 73, 282.
"""
function epb(day1::Real, day2::Real)
    1900.0 + ((day1 - JULIANDAY2000) + (day2 + BD1900))/DAYINYEAR1900
end

"""
    epb2jd(epoch::Real)

Besselian Epoch to Julian Date.

# Input

 - `epoch` -- Besselian Epoch (e.g. 1957.3)

# Output

 - `mjd0`  -- MJD zero-point: always 2400000.5
 - `mjd`   -- Modified Julian Date

# Note

The Julian Date is returned in two pieces, in the usual ERFA manner,
which is designed to preserve time resolution.  The Julian Date is
available as a single number by adding MJD0 and MJD.

# References

Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
"""
function epb2jd(epoch::Real)
    (mjd0 = MJDAY0, mjd = 15019.81352 + (epoch - 1900.0) * DAYINYEAR1900)
end

"""
    epj(day1::Real, day2::Real)

Julian Date to Julian Epoch.

# Input

 - `day1`  -- Julian Date (see note)
 - `day2`  -- Julian Date (see note)

# Output

 - `epoch` -- Julian Epoch

# Note

The Julian Date is supplied in two pieces, in the usual ERFA manner,
which is designed to preserve time resolution.  The Julian Date is
available as a single number by adding day1 and day2.  The maximum
resolution is achieved if day1 is 2451545.0 (J2000.0).

# References

Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
"""
function epj(day1::Real, day2::Real)
    2000.0 + ((day1 - JULIANDAY2000) + day2) / DAYINYEAR2000
end

"""
    epj2jd(epoch::Real)

Julian Epoch to Julian Date.

# Input

 - `epoch` -- Julian Epoch (e.g. 1996.8)

# Output

 - `mjd0`  -- MJD zero-point: always 2400000.5
 - `mjd`   -- Modified Julian Date

# Note

The Julian Date is returned in two pieces, in the usual ERFA manner,
which is designed to preserve time resolution.  The Julian Date is
available as a single number by adding MJD0 and MJD.

# References

Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
"""
function epj2jd(epoch::Real)
    (mjd0 = MJDAY0, mjd = MODJULDAY0 + (epoch - 2000.0) * DAYINYEAR2000)
end

"""
    jd2cal(day1::Real, day2::Real)

Julian Date to Gregorian year, month, day, and fraction of a day.

# Input

 - `day1`  -- Julian Date (Notes 1, 2)
 - `day2`  -- Julian Date (Notes 1, 2)

# Output

 - `year`  -- year
 - `month` -- month
 - `day`   -- day
 - `fraction` -- fraction of day

# Note

1) The earliest valid date is -68569.5 (-4900 March 1).  The largest
   value accepted is 1e9.

2) The Julian Date is apportioned in any convenient way between the
   arguments dj1 and dj2.  For example, JD=2450123.7 could be
   expressed in any of these ways, among others:

          day1             day2

        2450123.7           0.0       (JD method)
        2451545.0       -1421.3       (J2000 method)
        2400000.5       50123.2       (MJD method)
        2450123.5           0.2       (date & time method)

   Separating integer and fraction uses the "compensated summation"
   algorithm of Kahan-Neumaier to preserve as much precision as
   possible irrespective of the jd1+jd2 apportionment.

3) In early eras the conversion is from the "proleptic Gregorian
   calendar"; no account is taken of the date(s) of adoption of the
   Gregorian calendar, nor is the AD/BC numbering convention observed.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 12.92
(p604).

Klein, A., A Generalized Kahan-Babuska-Summation-Algorithm.
Computing, 76, 279-293 (2006), Section 3.
"""
function jd2cal(day1::Real, day2::Real)

    @assert JDMIN <= (day1 + day2) <= JDMAX "Day is out of range."

    # Separate day and fraction where fraction in range [-0.5, 0.5].
    day::Integer = convert(Int, round(day1)) + convert(Int, round(day2))
    fracs::Vector{AbstractFloat} = [(day1 - round(day1)), (day2 - round(day2))]
    
    # Compute frac1 + frac2 + 0.5 using compensated summation (Klein 2006).
    cs::Float64, s::Float64, t::Float64 = 0.0, 0.5, 0.0
    for x in fracs
        t = s + x
        cs += abs(s) >= abs(x) ? (s - t) + x : (x - t) + s
        s = t
        if s >= 1.0
            day += 1
            s -= 1.0
        end
    end
    frac::Float64 = s + cs
    cs = frac - s

    # Correct for negative fraction
    if frac < 0.0
        # Compensated summation: assume that |s| <= 1.0
        frac = s + 1.0
        cs += (1.0 - frac) + s
        s = frac
        frac = s + cs
        cs = frac - s
        day -= 1
    end

    # Correct for fraction that is >= 1.0 when rounded to float
    if frac - 1.0 >= -typemin(typeof(frac))/4.0
        # Compensated summation: assume that |s| <= 1.0
        t = s - 1.0
        cs += (s - t) - 1.0
        s = t
        frac = s + cs
        if -typemin(typeof(frac))/2.0 < frac
            day += 1
            frac = maximum((frac, 0.0))
        end
    end

#    day2cal(day)
    ll::Integer = day + 68569
    nn::Integer = (4 * ll)÷146097
    ll -= (146097 * nn + 3)÷4
    ii = (4000 * (ll + 1))÷1461001
    ll -= (1461 * ii)÷4 - 31
    kk::Integer = (80 * ll)÷2447
    day = ll - (2447 * kk)÷80
    ll = kk÷11
    month::Integer = kk + 2 - 12*ll
    year::Integer = 100 * (nn - 49) + ii + ll

    (year = year, month = month, day = day, fraction = frac)
end

"""
    jdcalf(ndp::Integer, day1::Real, day2::Real)

Julian Date to Gregorian Calendar, expressed in a form convenient for
formatting messages: rounded to a specified precision.

# Input

 - `ndp`   -- number of decimal places of days in fraction
 - `day1`  -- day1+day2 = Julian Date (Note 1)
 - `day2`  -- day1+day2 = Julian Date (Note 1)

# Output

 - `year`  -- year in Gregorian calendar
 - `month` -- month in Gregorian calendar
 - `day`   -- day in Gregorian calendar
 - `fraction` -- fraction in Gregorian calendar

# Note

1) The Julian Date is apportioned in any convenient way between the
   arguments dj1 and dj2.  For example, JD=2450123.7 could be
   expressed in any of these ways, among others:

          day1            day2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

2) In early eras the conversion is from the "Proleptic Gregorian
   Calendar"; no account is taken of the date(s) of adoption of the
   Gregorian Calendar, nor is the AD/BC numbering convention observed.

3) See also the function jd2cal.

4) The number of decimal places ndp should be 4 or less if internal
   overflows are to be avoided on platforms which use 16-bit integers.
"""
function jdcalf(ndp::Integer, day1::Real, day2::Real)

    # Denominator of fraction (e.g., 100 for 2 decimal places)
    jj::Integer, denom::Float64 = 0 <= ndp <= 9 ? (0, 10.0^ndp) : (1, 1.0)

    d1, d2 = abs(day1) >= abs(day2) ? (day1, day2) : (day2, day1)

    # Re-align the midnight (without roundng error)
    d1 -= 0.5

    # Separate day and fraction (as precisely as possible)
    f1, f2 = d1 - round(d1), d2 - round(d2)
    djd = round(d1) + round(d2)
    dd = round(f1 + f2)
    ff = (f1 - dd) + f2
    if ff < 0.0
        ff, dd = ff + 1.0, dd - 1.0
    end
    djd += dd

    # Round the total fraction to the specified number of places
    rf = round(ff*denom)/denom

    # Re-align to noon
    djd += 0.5

    # Convert to Gregorian calendar
    year::Integer, month::Integer, day::Integer, fraction::Float64 = jd2cal(djd, rf)
    # fraction = round(ff * denom)::Integer

    (year = year, month = month, day = day, fraction = convert(Integer, round(ff * denom)))
end
