#### Vector-Matrix / Angle Operations

"""
Decompose radians into degrees, arcminutes, arcseconds, and fraction.

# Arguments
- `npd::Integer`: number of useful digits
- `angle::Float64`: angle in radians

# Returns
- `dms::NamedTuple{(:sign, :degree, :minute, :second, :fraction)}`:
   angle in sign, degrees, minutes, seconds, and fraction

Notes:

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the size
   of angle, the format of Float64 on the target platform, and the risk
   of overflowing dms[3].  On a typical platform, for angle up to
   2pi, the available floating-point precision might correspond to
   ndp=12.  However, the practical limit is typically ndp=9, set by
   the capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of angle may exceed 2pi.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 2pi and rounds up to 2pi and rounds up to 360
   degree, by testing for dms[0]=360 and setting dms[0-3] to zero.
"""
function a2af(ndp::Integer, angle::Float64)
    NamedTuple{(:sign, :degree, :minute, :second, :fraction)}
    (values(d2tf(ndp, angle*15/2/pi)))
end

"""
Decompose radians into hours, minutes, seconds, and fraction.

# Arguments
- `npd::Integer`: number of useful digits
- `angle::Float64`: angle in radians

# Returns
- `hms::NamedTuple{(:sign, :hour, :minute, :second, :fraction)}`:
   angle in sign, hour, minutes, seconds, and fraction

Notes:

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the size
   of angle, the format of Float64 on the target platform, and the risk
   of overflowing hms[3].  On a typical platform, for angle up to
   2pi, the available floating-point precision might correspond to
   ndp=12.  However, the practical limit is typically ndp=9, set by
   the capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of angle may exceed 2pi.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 2pi and rounds up to 2pi and rounds up to 24
   hours, by testing for hms[0]=24 and setting hms[0-3] to zero.
"""
function a2tf(ndp::Integer, angle::Float64)
    NamedTuple{(:sign, :hour, :minute, :second, :fraction)}(d2tf(ndp, angle/2/pi))
end

"""
    Convert degrees, arcminutes, arcseconds to radians.

# Arguments
- `sign::Character`: sign of arc
- `degree::Integer`: degrees of arc
- `minute::Integer`: minutes of arc
- `second::Integer`: seconds of arc

# Returns
- `angle::Float64`: angle in radians

Notes:

1)  The result is computed even if any of the range checks fail.

2)  Negative ideg, iamin and/or asec produce a warning status, but
    the absolute value is used in the conversion.

3)  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.
"""
function af2a(sign::Char, degree::Integer, minute::Integer, second::Float64)
    @assert 0   <= degree < 360  "degree out of range [0-359]."
    @assert 0   <= minute <  60  "minute out of range [0-59]."
    @assert 0.0 <= second < 60.0 "second out of range [0-60]."
    
    deg2rad(1/3600) * (sign == '-' ? -1.0 : 1.0) *
       (60.0 * (60.0 * abs(degree) + abs(minute)) + abs(second))
end

"""
    Normalize angle into the range 0 <= a < 2p.

# Arguments
- `angle::Float64`: angle in radians

# Returns
- `angle::Float64`: angle in radians in range 0-2pi
"""
function anp(angle::Float64)
    mod2pi(angle)
end

"""
    Normalize angle into the range -pi <= a < +pi

# Arguments
- `angle::Float64`: angle in radians

# Returns
- `angle::Float64`: angle in radians in range +/-pi
"""
function anpm(angle::Float64)
    rem2pi(angle, RoundNearest)
end

"""
    d2tf(ndp, day)

Decompose days to sign, hours, minutes, seconds, fraction.

# Arguments
- `npd::Integer`: number of usefule digits
- `day::Float64`: interval in days

# Returns
- `hms::NamedTuple{(:sign, :hour, :minute, :second, :fraction)}`:
hms in sign, hours, minutes, seconds, and fraction

Notes:

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the size
   of days, the format of Float64 on the target platform, and the risk
   of overflowing hms[3].  On a typical platform, for days up to
   1.0, the available floating-point precision might correspond to
   ndp=12.  However, the practical limit is typically ndp=9, set by
   the capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of days may exceed 1.0.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 1.0 and rounds up to 24 hours, by testing for
   hms[0]=24 and setting hms[0-3] to zero.
"""
function d2tf(ndp::Integer, day::Float64)

    a = SECPERDAY * abs(day)
    if ndp < 0
        rs = prod(n == 2 || n == 4 ? 6 : 10 for n=1:-ndp)
        a = rs * round(a/rs)
    else
        rs = 10^ndp
    end

    rh, rm = 3600.0 * rs, 60.0 * rs

    sn = day >= 0.0 ? '+' : '-'
    a  = round(rs * a)
    ah = convert(Integer, trunc(a/rh))
    a -= ah*rh
    am = convert(Integer, trunc(a/rm))
    a -= am*rm
    as = convert(Integer, trunc(a/rs))
    af = convert(Integer, a - as*rs)

    NamedTuple{(:sign, :hour, :minute, :second, :fraction)}((sn, ah, am, as, af))
end

"""
    Convert hours, minutes, seconds to radians.

# Arguments
- `sign::Char`: sign:  '-' = negative, otherwise positive
- `hour::Integer`: hours
- `minute::Integer`: minutes
- `second::Integer`: seconds

# Returns
- `angle::Float64`: angle in radians

Notes:

1)  The result is computed even if any of the range checks fail.

2)  Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

3)  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.
"""
function tf2a(sign::Char, hour::Integer, minute::Integer, second::Float64)
    @assert 0   <= hour   <   24  "hour out of range [0-23]."
    @assert 0   <= minute <   60  "minute out of range [0-59]."
    @assert 0.0 <= second < 60.0 "second out of range [0-60]."
    
    15*deg2rad(1/3600) * (sign == '-' > -1.0 : 1.0) *
        (60.0 * (60.0 * abs(hour) + abs(minute)) + abs(second))
end

"""
    Convert hours, minutes, seconds to days.

# Arguments
- `sign::Char`: sign:  '-' = negative, otherwise positive
- `hour::Integer`: hours
- `minute::Integer`: minutes
- `second::Integer`: seconds

# Returns
- `day::Float64`: interval in days

Notes:

1)  The result is computed even if any of the range checks fail.

2)  Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

3)  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.
"""
function tf2d(sign::Char, hour::Integer, minute::Integer, second::Float64)
    @assert 0   <= hour   <   24  "hour out of range [0-23]."
    @assert 0   <= minute <   60  "minute out of range [0-59]."
    @assert 0.0 <= second < 60.0 "second out of range [0-60]."
    
    (sign == '-' > -1.0 : 1.0) *
        (60.0 * (60.0 * abs(hour) + abs(minute)) + abs(second)) / SECPERDAY
end

#### Vector - Matrix / Build Rotations

"""
    Rotate an r-matrix about the x-axis.

# Arguments:
- `ϕ::Float64`: angle (radians)
- `r::Matrix{Float64}`: r-matrix

# Returns:
- `r::Matrix{Float64}`: r-matrix, rotated

Notes:

1) Calling this function with positive ϕ incorporates in the
   supplied r-matrix r an additional rotation, about the x-axis,
   anticlockwise as seen looking towards the origin from positive x.

2) The additional rotation can be represented by this matrix:

       (  1        0            0      )
       (                               )
       (  0   + cos(ϕ)   + sin(ϕ)  )
       (                               )
       (  0   - sin(ϕ)   + cos(ϕ)  )

"""
function rx(ϕ::Float64, r::Matrix{Float64})
    Rx(ϕ)*r
end

"""
    Rotate an r-matrix about the y-axis.

# Arguments:
- `θ::Float64`: angle (radians)
- `r::Matrix{Float64}`: r-matrix

# Returns:
- `r::Matrix{Float64}`: r-matrix, rotated

Notes:

1) Calling this function with positive theta incorporates in the
   supplied r-matrix r an additional rotation, about the y-axis,
   anticlockwise as seen looking towards the origin from positive y.

2) The additional rotation can be represented by this matrix:

       (  + cos(θ)     0      - sin(θ)  )
       (                                        )
       (       0           1           0        )
       (                                        )
       (  + sin(θ)     0      + cos(θ)  )
"""
function ry(θ::Float64, r::Matrix{Float64})
    Ry(θ)*r
end

"""
    Rotate an r-matrix about the z-axis.

# Arguments:
- `ψ::Float64`: angle (radians)
- `r::Matrix{Float64}`: r-matrix

# Returns:
- `r::Matrix{Float64}`: r-matrix, rotated

Notes:

1) Calling this function with positive ψ incorporates in the
   supplied r-matrix r an additional rotation, about the z-axis,
   anticlockwise as seen looking towards the origin from positive z.

2) The additional rotation can be represented by this matrix:

       (  + cos(ψ)   + sin(ψ)     0  )
       (                              )
       (  - sin(ψ)   + cos(ψ)     0  )
       (                              )
       (       0            0      1  )

"""
function rz(ψ::Float64, r::Matrix{Float64})
    Rz(ψ)*r
end

#### Vector - Matrix / Copy, Extend, Extract

"""
    Copy a p-vector.

# Arguments:
- `p::Vector{Float64}`: p-vector to be copied

# Returns:
- `c::Vector{Float64}`: copy
"""
cp(p::Vector{Float64}) = copy(p)

"""
    Copy a position/velocity vector.

# Arguments:
- `pv::Vector{Vector{Float64}}`: position/velocity vector to be copied

# Returns:
- `c::Vector{Vector{Float64}}`: copy
"""
cpv(pv::Vector{Vector{Float64}}) = deepcopy(pv)

"""
    Copy an r-matrix.

# Arguments:
- `r::Matrix{Float64}`: r-matrix to be copied

# Returns:
- `c::Matrix{Float64}`: copy
"""
cr(r::Matrix{Float64}) = copy(r)

"""
    Extend a p-vector to a pv-vector by appending a zero velocity.

# Arguments:
- `p::Vector{Float64}`: p-vector

# Returns:
- `pv::Vector{Vector{Float64}}`: pv-vector
"""
p2pv(p::Vector{Float64}) = [p, zeros(Float64, 3)]

"""
    Discard velocity component of a pv-vector.

# Arguments:
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `p::Vector{Float64}`: p-vector
"""
pv2p(pv::Vector{Vector{Float64}}) = pv[1]

#### Vector - Matrix / Initialization

"""
Initialize an r-matrix to the identity matrix.

# Returns:
   r       Matrix{Float64}    r-matrix
"""
ir() = I + zeros(Float64, 3, 3)

"""
Zero a p-vector.

# Returns:
   p        Vector{Float64}      zero p-vector
"""
zp() = zeros(Float64, 3)

"""
Zero a pv-vector.

# Returns:
   pv       Vector{Vector{Float64}}      zero pv-vector
"""
zpv() = [zeros(Float64, 3), zeros(Float64, 3)]

"""
Initialize an r-matrix to the null matrix.

# Returns:
   r        Matrix{Float64}    r-matrix
"""
zr() = zeros(Float64, 3, 3)

#### Vector - Matrix / Matrix Operations

"""
    Multiply two r-matrices.

# Arguments
- `a::Matrix{Float64}`: first r-matrix
- `b::Matrix{Float64}`: second r-matrix

# Returns:
- `atb::Matrix{Float64}`: a * b

Note:

   It is permissible to re-use the same array for any of the
   arguments.
"""
rxr(a::Matrix{Float64}, b::Matrix{Float64}) = a*b

"""
    Transpose an r-matrix.

# Arguments
- `r::Matrix{Float64}`: r-matrix

# Returns:
- `rt::Matrix{Float64}`: transpose

Note:

   It is permissible for r and rt to be the same array.
"""
tr(r::Matrix{Float64}) = r'

#### Vector - Matrix / Matrix-Vector Products

"""
    Multiply a p-vector by an r-matrix.

# Arguments
- `r::Matrix{Float64}`: r-matrix
- `p::Vector{Float64}`: p-vector

# Returns:
- `rp::Vector{Float64}`: r * p

Note:
   It is permissible for p and rp to be the same array.
"""
rxp(r::Matrix{Float64}, p::Vector{Float64}) = r*p

"""
    Multiply a pv-vector by an r-matrix.

# Arguments
- `r::Matrix{Float64}`: r-matrix
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `rpv::Vector{Vector{Float64}}`: r * pv

Notes:

1) The algorithm is for the simple case where the r-matrix r is not a
   function of time.  The case where r is a function of time leads to
   an additional velocity component equal to the product of the
   derivative of r and the position vector.

2) It is permissible for pv and rpv to be the same array.
"""
rxpv(r::Matrix{Float64}, pv::Vector{Vector{Float64}}) = [r*pv[1], r*pv[2]]

"""
    Multiply a p-vector by the transpose of an r-matrix.

# Arguments
- `r::Matrix{Float64}`: r-matrix
- `p::Vector{Float64}`: p-vector

# Returns:
- `trp::Vector{Float64}`: r^T * p

Note:
   It is permissible for p and trp to be the same array.
"""
trxp(r::Matrix{Float64}, p::Vector{Float64}) = r'*p

"""
    Multiply a pv-vector by the transpose of an r-matrix.

# Arguments
- `r::Matrix{Float64}`: r-matrix
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `trpv::Vector{Vector{Float64}}`: r^T * pv

Notes:

1) The algorithm is for the simple case where the r-matrix r is not a
   function of time.  The case where r is a function of time leads to
   an additional velocity component equal to the product of the
   derivative of the transpose of r and the position vector.

2) It is permissible for pv and rpv to be the same array.
"""
trxpv(r::Matrix{Float64}, pv::Vector{Vector{Float64}}) = [r'*pv[1], r'*pv[2]]

#### Vector - Matrix / Rotation Vectors

"""
    Express an r-matrix as an r-vector.

# Arguments:
- `r::Matrix{Float64}`: rotation matrix

# Returns:
- `w::Vector{Float64}`: rotation vector (Note 1)

Notes:

1) A rotation matrix describes a rotation through some angle about
   some arbitrary axis called the Euler axis.  The "rotation vector"
   returned by this function has the same direction as the Euler axis,
   and its magnitude is the angle in radians.  (The magnitude and
   direction can be separated by means of the function eraPn.)

2) If r is null, so is the result.  If r is not a rotation matrix the
   result is undefined; r must be proper (i.e. have a positive
   determinant) and real orthogonal (inverse = transpose).

3) The reference frame rotates clockwise as seen looking along the
   rotation vector from the origin.
"""
function rm2v(r::Matrix{Float64})
    x, y, z = r[2,3] - r[3,2], r[3,1] - r[1,3], r[1,2] - r[2,1]
    s2, c2 = norm([x, y, z]), r[1,1] + r[2,2] + r[3,3] - 1
    s2 > 0 ? [x, y, z]*atan(s2, c2)/s2 : [0.0, 0.0, 0.0]
end

"""
    Form the r-matrix corresponding to a given r-vector.

# Arguments:
- `w::Vector{Float64}`: rotation vector (Note 1)

# Returns:
- `r::Matrix{Float64}`: rotation matrix

Notes:

1) A rotation matrix describes a rotation through some angle about
   some arbitrary axis called the Euler axis.  The "rotation vector"
   supplied to This function has the same direction as the Euler axis,
   and its magnitude is the angle in radians.

2) If w is null, the identity matrix is returned.

3) The reference frame rotates clockwise as seen looking along the
   rotation vector from the origin.
"""

function rv2m(w::Vector{Float64})
    #  Euler angle (magnitude of rotation vector)
    ϕ = norm(w)
    #  Euler axis (direction of rotation vector), perhaps null
    k = -(ϕ > 0 ? w/ϕ : w)
    I + sin(ϕ)*vec2mat(k) + (1-cos(ϕ))*vec2mat(k)*vec2mat(k)
end

#### Vector - Matrix / Separation and Angle

"""
    Position-angle from two p-vectors.

# Arguments:
- `a::Vector{Float64}`: direction of reference point
- `b::Vector{Float64}`: direction of point whose PA is required

# Returns:
- `θ::Float64`: position angle of b with respect to a (radians)

Notes:

1) The result is the position angle, in radians, of direction b with
   respect to direction a.  It is in the range -pi to +pi.  The sense
   is such that if b is a small distance "north" of a the position
   angle is approximately zero, and if b is a small distance "east" of
   a the position angle is approximately +pi/2.

2) The vectors a and b need not be of unit length.

3) Zero is returned if the two directions are the same or if either
   vector is null.

4) If vector a is at a pole, the result is ill-defined.
"""
function pap(a::Vector{Float64}, b::Vector{Float64})
    if norm(a) == 0 || norm(b) == 0
        θ = 0.0
    else
        #  The north axis tangent from a (arbitrary length)
        η = [-a[1]*a[3], -a[2]*a[3], sum(a[1:2].^2)]
        #  The east axis tanget from a (same length)
        ξ = vec2mat(η)*a/norm(a)
        # Resolve into components along the north and east axes
        θ = (b.-a)'*ξ == 0 && (b.-a)'*η == 0 ? 0.0 : atan((b.-a)'*ξ, (b.-a)'*η)
    end
    θ
end

"""
    Position-angle from spherical coordinates.

# Arguments:
- `λa::Float64`: longitude of point A (e.g. RA) in radians
- `ϕa::Float64`: latitude of point A (e.g. Dec) in radians
- `λb::Float64`: longitude of point B
- `ϕb::Float64`: latitude of point B

# Returns:
- `θ::Float64`: position angle of B with respect to A

Notes:

1) The result is the bearing (position angle), in radians, of point B
   with respect to point A.  It is in the range -pi to +pi.  The sense
   is such that if B is a small distance "east" of point A, the
   bearing is approximately +pi/2.

2) Zero is returned if the two points are coincident.
"""
function pas(λa::Float64, ϕa::Float64, λb::Float64, ϕb::Float64)
    x = sin(ϕb)*cos(ϕa) - cos(ϕb)*sin(ϕa)*cos(λb - λa)
    y = sin(λb - λa)*cos(ϕb)
    x != 0 || y != 0 ? atan(y, x) : 0.0
end

"""
    Angular separation between two p-vectors.

# Arguments:
- `a::Vector{Float64}`: first p-vector (not necessarily unit length)
- `b::Vector{Float64}`: second p-vector (not necessarily unit length)

# Returns:
- `θ::Float64`: angular separation (radians, always positive)

Notes:

1) If either vector is null, a zero result is returned.

2) The angular separation is most simply formulated in terms of scalar
   product.  However, this gives poor accuracy for angles near zero
   and pi.  The present algorithm uses both cross product and dot
   product, to deliver full accuracy whatever the size of the angle.
"""
function sepp(a::Vector{Float64}, b::Vector{Float64})
    #  Sine of angle between the vectors, multiplied by the two moduli
    #  Cosine of the angle, multiplied by the two moduli
    cosθ, sinθ = sum(a.*b), norm(vec2mat(a)*b)
    sinθ != 0 || cosθ != 0 ? atan(sinθ, cosθ) : 0.0
end

"""
    Angular separation between two sets of spherical coordinates.

# Arguments:
- `λa::Float64`: first longitude (radians)
- `ϕa::Float64`: first latitude (radians)
- `λb::Float64`: second longitude (radians)
- `ϕb::Float64`: second latitude (radians)

# Returns:
- `θ::Float64`: angular separation (radians)

"""
function seps(λa::Float64, ϕa::Float64, λb::Float64, ϕb::Float64)
    #=
    #  Spherical to Cartesian
    a = [cos(λa)*cos(ϕa), sin(λa)*cos(ϕa), sin(ϕa)]
    b = [cos(λb)*cos(ϕb), sin(λb)*cos(ϕb), sin(ϕb)]
    #  Sine of angle between the vectors, multiplied by the two moduli
    sinθ = norm([0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]*b)
    #  Cosine of the angle, multiplied by the two moduli
    cosθ = sum(a.*b)
    sinθ != 0 || cosθ != 0 ? atan(sinθ, cosθ) : 0.0
    =#
    @inline sepp(s2c(λa, ϕa), s2c(λb, ϕb))
end

#### Vector - Matrix / Spherical-Cartesian

"""
    P-vector to spherical coordinates.

# Arguments:
- `p::Vector{Float64}`: p-vector

# Returns:
- `θ::Float64`: longitude angle (radians)
- `ϕ::Float64`: latitude angle (radians)

Notes:

1) The vector p can have any magnitude; only its direction is used.

2) If p is null, zero θ and ϕ are returned.

3) At either pole, zero θ is returned.
"""
function c2s(pos::Vector{Float64})
    NamedTuple{(:θ, :ϕ)}
    ((sum(pos[1:2].^2) == 0 ? 0 : atan(pos[2], pos[1]),
      pos[3] == 0 ? 0 : atan(pos[3], norm(pos[1:2]))))
end

"""
    P-vector to spherical polar coordinates.

# Arguments:
- `p::Vector{Float64}`: p-vector

# Returns:
- `θ::Float64`: longitude angle (radians)
- `ϕ::Float64`: latitude angle (radians)
- `r::Float64`: radial distance

Notes:

1) If P is null, zero θ, ϕ and r are returned.

2) At either pole, zero θ is returned.

"""
function p2s(pos::Vector{Float64})
    @inline NamedTuple{(:θ, :ϕ, :r)}((c2s(pos)..., norm(pos)))
end

"""
    Convert position/velocity from Cartesian to spherical coordinates.

# Arguments:
- `posvel::Vector{Vector{Float64}}`: position-velocity-vector

# Returns:
- `θ::Float64`: longitude angle (radians)
- `ϕ::Float64`: latitude angle (radians)
- `r::Float64`: radial distance
- `dθ::Float64`: rate of change of θ
- `dϕ::Float64`: rate of change of ϕ
- `dr::Float64`: rate of change of r

Notes:

1) If the position part of pv is null, theta, ϕ, td and pd are
   indeterminate.  This is handled by extrapolating the position
   through unit time by using the velocity part of pv.  This moves the
   origin without changing the direction of the velocity component.
   If the position and velocity components of pv are both null, zeroes
   are returned for all six results.

2) If the position is a pole, theta, td and pd are indeterminate.  In
   such cases zeroes are returned for all three.
"""
function pv2s(pv::Vector{Vector{Float64}})
    
    x,  y,  z  = pv[(sum(pv[1].^2) == 0) ? 2 : 1]
    dx, dy, dz = pv[2]

    if sum(pv[1][1:2].^2) != 0
        θ  = atan(y, x)
        ϕ  = atan(z, norm([x, y]))
        dθ = (x*dy - y*dx) / (x*x+y*y)
        dϕ = (dz*(x*x+y*y) - z*(x*dx+y*dy)) / (sum([x, y, z].^2)*norm([x, y]))
    else
        θ,  ϕ  = 0, (z != 0) ? atan(z, norm([x, y])) : 0
        dθ, dϕ = 0, 0
    end
    r  = norm(pv[1])
    dr = norm([x, y, z]) != 0 ? (x*dx+y*dy+z*dz)/norm([x, y, z]) : 0
    NamedTuple{(:θ, :ϕ, :r, :dθ, :dϕ, :dr)}((θ, ϕ, r, dθ, dϕ, dr))
end

"""
    Convert spherical coordinates to Cartesian.

# Arguments:
- `θ::Float64`: longitude angle (radians)
- `ϕ::Float64`: latitude angle (radians)

# Returns:
- `c::Vector{Float64}`: direction cosines

"""
function s2c(θ::Float64, ϕ::Float64)
    [cos(θ)*cos(ϕ), sin(θ)*cos(ϕ), sin(ϕ)]
end

"""
    Convert spherical polar coordinates to p-vector.

# Arguments:
- `θ::Float64`: longitude angle (radians)
- `ϕ::Float64`: latitude angle (radians)
- `r::Float64`: radial distance

# Returns:
- `p::Vector{Float64}`: Cartesian coordinates
"""
s2p(θ::Float64, ϕ::Float64, r::Float64) = r*s2c(θ, ϕ)

"""
    Convert position/velocity from spherical to Cartesian coordinates.

# Arguments:
- `θ::Float64`: longitude angle (radians)
- `ϕ::Float64`: latitude angle (radians)
- `r::Float64`: radial distance
- `dθ::Float64`: rate of change of θ
- `dϕ::Float64`: rate of change of ϕ
- `dr::Float64`: rate of change of r

# Returns:
- `pv::Vector{Vector{Float64}}`: pv-vector
"""
function s2pv(θ::Float64, ϕ::Float64, r::Float64,
              dθ::Float64, dϕ::Float64, dr::Float64)
    [[r*cos(θ)*cos(ϕ), r*sin(θ)*cos(ϕ), r*sin(ϕ)],
     [-r*dθ*sin(θ)*cos(ϕ) - cos(θ)*(r*dϕ*sin(ϕ) - dr*cos(ϕ)),
      r*dθ*cos(θ)*cos(ϕ) - sin(θ)*(r*dϕ*sin(ϕ) - dr*cos(ϕ)),
      r*dϕ*cos(ϕ) + dr*sin(ϕ)]]
end

#### Vector - Matrix / Vector Operations

"""
    P-vector inner (=scalar=dot) product.

# Arguments:
- `a::Vector{Float64}`: first p-vector
- `b::Vector{Float64}`: second p-vector

# Returns:
- `r::Float64        a . b
"""
pdp(a::Vector{Float64}, b::Vector{Float64}) = sum(a.*b)

"""
    Modulus of p-vector.

# Arguments:
- `p::Vector{Float64}`: p-vector

# Returns:
- `r::Float64`: modulus
"""
pm(p::Vector{Float64}) = norm(p)

"""
    P-vector subtraction.

# Arguments:
- `a::Vector{Float64}`: first p-vector
- `b::Vector{Float64}`: second p-vector

# Returns:
- `amb::Vector{Float64}`: a - b
"""
pmp(a::Vector{Float64}, b::Vector{Float64}) = a .- b

"""
    Convert a p-vector into modulus and unit vector.

# Arguments:
- `p::Vector{Float64}`: p-vector

# Returns:
- `r::Float64`: modulus
- `u::Vector{Float64}`: unit vector

Notes:

1) If p is null, the result is null.  Otherwise the result is a unit
   vector.
"""
function pn(p::Vector{Float64})
    NamedTuple{(:modulus, :unit)}
    (norm(p) == 0 ? (0., [0., 0., 0.]) : (norm(p), p./norm(p)))
end

"""
    P-vector addition.

# Arguments:
- `a::Vector{Float64}`: first p-vector
- `b::Vector{Float64}`: second p-vector

# Returns:
- `apb::Vector{Float64}`: a + b
"""
ppp(a::Vector{Float64}, b::Vector{Float64}) = a.+b

"""
    P-vector plus scaled p-vector.

# Arguments:
- `a::Vector{Float64}`: first p-vector
- `s::Float64`: scalar (multiplier for b)
- `b::Vector{Float64}`: second p-vector

# Returns:
- `apsb::Vector{Vector{Float64}}`: a + s*b
"""
ppsp(a::Vector{Float64}, s::Float64, b::Vector{Float64}) = a + s*b

"""
    Inner (=scalar=dot) product of two pv-vectors.

# Arguments:
- `a::Matrix{Float64}`: first pv-vector
- `b::Matrix{Float64}`: second pv-vector

# Returns:
- `adb::Vector{Float64}`: a . b (see note)

Note:

   If the position and velocity components of the two pv-vectors are (
   ap, av ) and ( bp, bv ), the result, a . b, is the pair of numbers
   ( ap . bp , ap . bv + av . bp ).  The two numbers are the
   dot-product of the two p-vectors and its derivative.
"""
function pvdpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})
   [sum(a[1].*b[1]), sum(a[1].*b[2] .+ a[2].*b[1])]
end

"""
    Modulus of pv-vector.

# Arguments:
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `r::Float64`: modulus of position component
- `s::Float64`: modulus of velocity component
"""
pvm(pv::Vector{Vector{Float64}}) = sqrt.(sum.([pv[1].^2, pv[2].^2]))

"""
    Subtract one pv-vector from another.

# Arguments:
- `a::Vector{Vector{Float64}}`: first pv-vector
- `b::Vector{Vector{Float64}}`: second pv-vector

# Returns:
- `amb::Vector{Vector{Float64}}`: a - b
"""
pvmpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}}) = [a[1] .- b[1], a[2] .- b[2]]

"""
    Add one pv-vector to another.

# Arguments:
- `a::Vector{Vector{Float64}}`: first pv-vector
- `b::Vector{Vector{Float64}}`: second pv-vector

# Returns:
- `apb::Vector{Vector{Float64}}`: a + b
"""
pvppv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}}) = [a[1] .+ b[1], a[2] .+ b[2]]

"""
    Update a pv-vector.

# Arguments:
- `dt::Float64`: time interval
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `upv::Vector{Vector{Float64}}`: p updated, v unchanged

Notes:

1) "Update" means "refer the position component of the vector to a new
   date dt time units from the existing date".

2) The time units of dt must match those of the velocity.
"""
pvu(dt::Float64, pv::Vector{Vector{Float64}}) = [pv[1] .+ dt.*pv[2], pv[2]]

"""
    Update a pv-vector, discarding the velocity component.

# Arguments:
- `dt::Float64`: time interval
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `p::Vector{Float64}`: p-vector

Notes:

1) "Update" means "refer the position component of the vector to a new
   date dt time units from the existing date".

2) The time units of dt must match those of the velocity.
"""
pvup(dt::Float64, pv::Vector{Vector{Float64}}) = pv[1] .+ dt*pv[2]

"""
    Outer (=vector=cross) product of two pv-vectors.

# Arguments:
- `a::Vector{Vector{Float64}}`: first pv-vector
- `b::Vector{Vector{Float64}}`: second pv-vector

# Returns:
- `axb::Matrix{Float64`: a x b

Notes:

1) If the position and velocity components of the two pv-vectors are (
   ap, av ) and ( bp, bv ), the result, a x b, is the pair of vectors
   ( ap x bp, ap x bv + av x bp ).  The two vectors are the
   cross-product of the two p-vectors and its derivative.
"""
function pvxpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})
    [vec2mat(a[1])*b[1], vec2mat(a[1])*b[2] .+ vec2mat(a[2])*b[1]]
end

"""
    P-vector outer (=vector=cross) product.

# Arguments:
- `a::Vector{Float64}`: first p-vector
- `b::Vector{Float64}`: second p-vector

# Returns:
- `axb::Vector{Float64}`: a x b
"""
pxp(a::Vector{Float64}, b::Vector{Float64}) = vec2mat(a)*b

"""
    Multiply a pv-vector by two scalars.

# Arguments:
- `s1::Float64`: scalar to multiply position component by
- `s2::Float64`: scalar to multiply velocity component by
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `spv::Vector{Vector{Float64}}`: pv-vector: p scaled by s1, v scaled by s2
"""
s2xpv(s1::Float64, s2::Float64, pv::Vector{Vector{Float64}}) = [s1*pv[1], s2*pv[2]]

"""
    Multiply a p-vector by a scalar.

# Arguments:
- `s::Float64`: scalar
- `p::Vector{Float64}`: p-vector

# Returns:
- `sp::Vector{Float64}`: s * p
"""
sxp(s::Float64, p::Vector{Float64}) = s*p

"""
    Multiply a pv-vector by a scalar.

# Arguments:
- `s::Float64`: scalar
- `pv::Vector{Vector{Float64}}`: pv-vector

# Returns:
- `spv::Vector{Vector{Float64}}`: s * pv
"""
sxpv(s::Float64, pv::Vector{Vector{Float64}}) = [s*pv[1], s*pv[2]]
