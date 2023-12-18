#### Astronomy / Geodetic & Geocentric Coordinates

"""
    eform(model::Symbol)

Earth reference ellipsoids.

# Input

 - `model` -- ellipsoid identifier (Note 1)

# Output

 - `radius` -- equatorial radius (meters, Note 2)
 - `oblate` -- oblateness (Note 2)

# Note

1) The identifier n is a number that specifies the choice of reference
   ellipsoid.  The following are supported:

      ellipsoid

      :WGS84
      :GRS80
      :WGS72

   The n value has no significance outside the ERFA software.  For
   convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.

2) The ellipsoid parameters are returned in the form of equatorial
   radius in meters (a) and flattening (f).  The latter is a number
   around 0.00335, i.e. around 1/298.

3) For the case where an unsupported n value is supplied, zero a and f
   are returned, as well as error status.

# References

Department of Defense World Geodetic System 1984, National Imagery and
Mapping Agency Technical Report 8350.2, Third Edition, p3-2.

Moritz, H., Bull. Geodesique 66-2, 187 (1992).

The Department of Defense World Geodetic System 1972, World Geodetic
System Committee, May 1974.

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), p220.
"""
function eform(model::Symbol)
    @assert model in [:WGS84, :GRS80, :WGS72] "Model not one of (:WGS72, :GRS80, :WGS84)."
    if model == :WGS84
        radius, oblate = wgs84_radius, wgs84_oblate
    elseif model == :GRS80
        radius, oblate = grs80_radius, grs80_oblate
    elseif model == :WGS72
        radius, oblate = wgs72_radius, wgs72_oblate
    end
    NamedTuple{(:radius, :oblate)}((radius, oblate))
end

"""
    gc2gd(model::Symbol, pos::Vector{Float64})
    
Transform geocentric coordinates to geodetic using the specified
reference ellipsoid.

# Input

 - `model::Int`: ellipsoid identifier (Note 1)
 - `pos::Vector{Float64}`: geocentric vector (Note 2)

# Output

 - `ϵ::Float64`: longitude (radians, east +ve, Note 3)
 - `ϕ::Float64`: latitude (geodetic, radians, Note 3)
 - `r::Float64`: height above ellipsoid (geodetic, Notes 2,3)

# Note

1) The identifier n is a number that specifies the choice of reference
   ellipsoid.  The following are supported:

      ellipsoid

      :WGS84
      :GRS80
      :WGS72

   The n value has no significance outside the ERFA software.  For
   convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.

2) The geocentric vector (xyz, given) and height (height, returned)
   are in meters.

3) An error status -1 means that the identifier n is illegal.  An
   error status -2 is theoretically impossible.  In all error cases,
   all three results are set to -1e9.

4) The inverse transformation is performed in the function eraGd2gc.
"""
function gc2gd(model::Symbol, pos::Vector{Float64})
    gc2gde(values(eform(model))..., pos)
end

"""
    gc2gde(radius::Float64, oblate::Float64, pos::Vector{Float64})

Transform geocentric coordinates to geodetic for a reference
ellipsoid of specified form.

# Input

 - `radius` -- equatorial radius (Notes 2,4)
 - `oblate` -- oblateness (Note 3)
 - `pos`    -- geocentric vector (Note 4)

# Output

 - `ϵ`     -- longitude (radians, east +ve)
 - `ϕ`     -- latitude (geodetic, radians)
 - `r`     -- height above ellipsoid (geodetic, Note 4)

# Note

1) This function is based on the GCONV2H Fortran subroutine by Toshio
   Fukushima (see reference).

2) The equatorial radius, a, can be in any units, but meters is the
   conventional choice.

3) The flattening, f, is (for the Earth) a value around 0.00335,
   i.e. around 1/298.

4) The equatorial radius, a, and the geocentric vector, xyz, must be
   given in the same units, and determine the units of the returned
   height, height.

5) If an error occurs (status < 0), elong, phi and height are
   unchanged.

6) The inverse transformation is performed in the function eraGd2gce.

7) The transformation for a standard ellipsoid (such as ERFA_WGS84)
   can more conveniently be performed by calling eraGc2gd, which uses
   a numerical code to identify the required A and F values.

# References

Fukushima, T., "Transformation from Cartesian to geodetic coordinates
accelerated by Halley's method", J.Geodesy (2006) 79: 689-693
"""
function gc2gde(radius::Float64, oblate::Float64, pos::Vector{Float64})
    @assert 0.0 <= oblate < 1.0 "Oblateness out of range [0 - 1)."
    @assert radius > 0.0 "Radius is <= 0."
    @assert (1.0 - (2.0 - oblate)*oblate) >= 0.0 "Oblateness is too larage."

    if sum(pos[1:2].^2) > 1e-32*radius^2
        #  Prepare Newton corection factors
        s ,  e = abs(pos[3])/radius, (2.0 - oblate)*oblate
        ec, pn = sqrt(1.0 - e), sqrt(sum(pos[1:2].^2))/radius
        a0 = sqrt((ec*pn)^2 + s^2)
        d0, f0 = ec*s*a0^3 + e*s^3, pn*a0^3 - e*(ec*pn)^3

        #  Prepare Halley correction factors
        b0 = 1.5*(e*s*ec*pn)^2*pn*(a0 - ec)
        s1 = d0*f0 - b0*s
        cc = ec*(f0*f0 - b0*ec*pn)
        ϵ = atan(pos[2], pos[1])
        ϕ = atan(s1/cc)
        r = (sqrt(sum(pos[1:2].^2))*cc + abs(pos[3])*s1 -
             radius*sqrt((ec*s1)^2 + cc^2))/sqrt(s1^2 + cc^2)
    else
        ϵ, ϕ, r = 0.0, π/2, abs(pos[3]) - radius*ec
    end
    NamedTuple{(:ϵ, :ϕ, :r)}((ϵ, ϕ, r))
end

"""
    gd2gc(model::Symbol, ϵ::Float64, ϕ::Float64, r::Float64)

Transform geodetic coordinates to geocentric using the specified
reference ellipsoid.

# Input

 - `model` -- ellipsoid identifier (Note 1)
 - `ϵ`     -- longitude (radians, east +ve)
 - `ϕ`     -- latitude (geodetic, radians, Note 3)
 - `r`     -- height above ellipsoid (geodetic, Notes 2,3)

# Output

 - `pos`   -- geocentric vector (Note 2)

# Note

1) The identifier n is a number that specifies the choice of reference
   ellipsoid.  The following are supported:

      ellipsoid

      :WGS84
      :GRS80
      :WGS72

   The n value has no significance outside the ERFA software.  For
   convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.

2) The height (height, given) and the geocentric vector (xyz,
   returned) are in meters.

3) No validation is performed on the arguments elong, phi and height.
   An error status -1 means that the identifier n is illegal.  An
   error status -2 protects against cases that would lead to
   arithmetic exceptions.  In all error cases, xyz is set to zeros.

4) The inverse transformation is performed in the function eraGc2gd.
"""
function gd2gc(model::Symbol, ϵ::Float64, ϕ::Float64, r::Float64)
    @inline gd2gce(values(eform(model))..., ϵ, ϕ, r)
end

"""
    gd2gce(radius::Float64, oblate::Float64, ϵ::Float64, ϕ::Float64, r::Float64)

Transform geodetic coordinates to geocentric for a reference ellipsoid
of specified form.

# Input

 - `radius` -- equatorial radius (Notes 1,4)
 - `oblate` -- oblateness (Notes 2,4)
 - `ϵ`      -- longitude (radians, east +ve)
 - `ϕ`      -- latitude (geodetic, radians, Note 4)
 - `r`      -- height above ellipsoid (geodetic, Notes 3,4)

# Output

 - `pos`    -- geocentric vector (Note 3)

# Note

1) The equatorial radius, a, can be in any units, but meters is the
   conventional choice.

2) The flattening, f, is (for the Earth) a value around 0.00335,
   i.e. around 1/298.

3) The equatorial radius, a, and the height, height, must be given in
   the same units, and determine the units of the returned geocentric
   vector, xyz.

4) No validation is performed on individual arguments.  The error
   status -1 protects against (unrealistic) cases that would lead to
   arithmetic exceptions.  If an error occurs, xyz is unchanged.

5) The inverse transformation is performed in the function eraGc2gde.

6) The transformation for a standard ellipsoid (such as ERFA_WGS84)
   can more conveniently be performed by calling eraGd2gc, which uses
   a numerical code to identify the required a and f values.

# References

Green, R.M., Spherical Astronomy, Cambridge University Press, (1985)
Section 4.5, p96.

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 4.22, p202.
"""
function gd2gce(radius::Float64, oblate::Float64, ϵ::Float64, ϕ::Float64, r::Float64)
    @assert (cos(ϕ)^2 + ((1.0-oblate)*sin(ϕ))^2) >= 0.0 "Equatorial radius is <= 0"
    d  = sqrt(cos(ϕ)^2 + ((1.0-oblate)*sin(ϕ))^2)
    rr = radius/d + r
    [rr*cos(ϕ)*cos(ϵ), rr*cos(ϕ)*sin(ϵ), ((1.0-oblate)^2*radius/d + r)*sin(ϕ)]
end
