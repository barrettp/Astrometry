#### Vector-Matrix / Angle Operations

"""
    a2af(ndp::Integer, angle::Float64)

Decompose radians into degrees, arcminutes, arcseconds, and fraction.

# Input

 - `npd`   --number of useful digits
 - `angle` -- angle in radians

# Output

 - `dms`   -- angle in sign, degrees, minutes, seconds, and fraction

# Note

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
   of angle, the format of Float64 on the target platform, and the
   risk of overflowing dms[3].  On a typical platform, for angle up to
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
    a2tf(ndp::Integer, angle::Float64)

Decompose radians into hours, minutes, seconds, and fraction.

# Input

 - `npd`   -- number of useful digits
 - `angle` -- angle in radians

# Output
 - `hms`   -- angle in sign, hour, minutes, seconds, and fraction

# Note

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
   of angle, the format of Float64 on the target platform, and the
   risk of overflowing hms[3].  On a typical platform, for angle up to
   2pi, the available floating-point precision might correspond to
   ndp=12.  However, the practical limit is typically ndp=9, set by
   the capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of angle may exceed 2pi.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 2pi and rounds up to 2pi and rounds up to 24
   hours, by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
"""
function a2tf(ndp::Integer, angle::Float64)
    NamedTuple{(:sign, :hour, :minute, :second, :fraction)}(d2tf(ndp, angle/2/pi))
end

"""
    af2a(sign::Char, degree::Integer, minute::Integer, second::Float64)

Convert degrees, arcminutes, arcseconds to radians.

# Input

 - `sign`   -- sign of arc
 - `degree` -- degrees of arc
 - `minute` -- minutes of arc
 - `second` -- seconds of arc

# Output

 - `angle`  -- angle in radians

# Note

1) The result is computed even if any of the range checks fail.

2) Negative ideg, iamin and/or asec produce a warning status, but the
   absolute value is used in the conversion.

3) If there are multiple errors, the status value reflects only the
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
    anp(angle::Float64)

Normalize angle into the range 0 <= a < 2p.

# Input

 - `angle` -- angle in radians

# Output

 - `angle` -- angle in radians in range 0-2pi
"""
function anp(angle::Float64)
    mod2pi(angle)
end

"""
    anpm(angle::Float64)

Normalize angle into the range -pi <= a < +pi

# Input

 - `angle` -- angle in radians

# Output

 - `angle` -- angle in radians in range +/-pi
"""
function anpm(angle::Float64)
    rem2pi(angle, RoundNearest)
end

"""
    d2tf(ndp::Integer, day::Float64)

Decompose days to sign, hours, minutes, seconds, fraction.

# Input

 - `npd`   -- number of usefule digits
 - `day`   -- interval in days

# Output

 - `hms`   -- hms in sign, hours, minutes, seconds, and fraction

# Note

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
   of overflowing hms[3].  On a typical platform, for days up to 1.0,
   the available floating-point precision might correspond to ndp=12.
   However, the practical limit is typically ndp=9, set by the
   capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of days may exceed 1.0.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 1.0 and rounds up to 24 hours, by testing for
   ihmsf[0]=24 and setting ihmsf[0-3] to zero.
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
    tf2a(sign::Char, hour::Integer, minute::Integer, second::Float64)

Convert hours, minutes, seconds to radians.

# Input

 - `sign`   -- sign:  '-' = negative, otherwise positive
 - `hour`   -- hours
 - `minute` -- minutes
 - `second` -- seconds

# Output

 - `angle`  -- angle in radians

# Note

1) The result is computed even if any of the range checks fail.

2) Negative ihour, imin and/or sec produce a warning status, but the
   absolute value is used in the conversion.

3) If there are multiple errors, the status value reflects only the
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
    tf2d(sign::Char, hour::Integer, minute::Integer, second::Float64)

Convert hours, minutes, seconds to days.

# Input

 - `sign`   -- sign:  '-' = negative, otherwise positive
 - `hour`   -- hours
 - `minute` -- minutes
 - `second` -- seconds

# Output

 - `day`    -- interval in days

# Note

1) The result is computed even if any of the range checks fail.

2) Negative ihour, imin and/or sec produce a warning status, but the
   absolute value is used in the conversion.

3) If there are multiple errors, the status value reflects only the
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
    rx(ϕ::Float64, r::Matrix{Float64})

Rotate an r-matrix about the x-axis.

# Input

 - `ϕ`     -- angle (radians)
 - `r`     -- r-matrix

# Output

 - `r`     -- r-matrix, rotated

# Note

1) Calling this function with positive ϕ incorporates in the supplied
   r-matrix r an additional rotation, about the x-axis, anticlockwise
   as seen looking towards the origin from positive x.

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
    ry(θ::Float64, r::Matrix{Float64})

Rotate an r-matrix about the y-axis.

# Input

 - `θ`     -- angle (radians)
 - `r`     -- r-matrix

# Output

 - `r`     -- r-matrix, rotated

# Note

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
    rz(ψ::Float64, r::Matrix{Float64})

Rotate an r-matrix about the z-axis.

# Input

 - `ψ`     -- angle (radians)
 - `r `    -- r-matrix

# Output

 - `r`     -- r-matrix, rotated

# Note

1) Calling this function with positive ψ incorporates in the supplied
   r-matrix r an additional rotation, about the z-axis, anticlockwise
   as seen looking towards the origin from positive z.

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
    cp(p::Vector{Float64})

Copy a p-vector.

# Input

 - `p`     -- p-vector to be copied

# Output

 - `c`     -- copy
"""
cp(p::Vector{Float64}) = copy(p)

"""
    cpv(pv::Vector{Vector{Float64}})

Copy a position/velocity vector.

# Input

 - `pv`     -- pv-vector to be copied

# Output

 - `c`      -- copy
"""
cpv(pv::Vector{Vector{Float64}}) = deepcopy(pv)

"""
    cr(r::Matrix{Float64})

Copy an r-matrix.

# Input

 - `r`     -- r-matrix to be copied

# Output

 - `c`     -- copy
"""
cr(r::Matrix{Float64}) = copy(r)

"""
    p2pv(p::Vector{Float64})

Extend a p-vector to a pv-vector by appending a zero velocity.

# Input

 - `p`     -- p-vector

# Output

 - `pv`    -- pv-vector
"""
p2pv(p::Vector{Float64}) = [p, zeros(Float64, 3)]

"""
    pv2p(pv::Vector{Vector{Float64}})

Discard velocity component of a pv-vector.

# Input

 - `pv`    -- pv-vector

# Output

 - `p`     -- p-vector
"""
pv2p(pv::Vector{Vector{Float64}}) = pv[1]

#### Vector - Matrix / Initialization

"""
    ir()

Initialize an r-matrix to the identity matrix.

# Output

 - `r`     -- r-matrix
"""
ir() = I + zeros(Float64, 3, 3)

"""
    zp()

Zero a p-vector.

# Output
 - `p`     -- zero p-vector
"""
zp() = zeros(Float64, 3)

"""
    zpv()

Zero a pv-vector.

# Output

 - `pv`    -- zero pv-vector
"""
zpv() = [zeros(Float64, 3), zeros(Float64, 3)]

"""
    zr()

Initialize an r-matrix to the null matrix.

# Output
 - `r`     -- r-matrix
"""
zr() = zeros(Float64, 3, 3)

#### Vector - Matrix / Matrix Operations

"""
    rxr(a::Matrix{Float64}, b::Matrix{Float64})

Multiply two r-matrices.

# Input

 - `a`     -- first r-matrix
 - `b`     -- second r-matrix

# Output

 - `atb`   -- a * b

# Note

1) It is permissible to re-use the same array for any of the
   arguments.
"""
rxr(a::Matrix{Float64}, b::Matrix{Float64}) = a*b

"""
    tr(r::Matrix{Float64})

Transpose an r-matrix.

# Input

 - `r`     -- r-matrix

# Output

 - `rt`    -- transpose

# Note

1) It is permissible for r and rt to be the same array.
"""
tr(r::Matrix{Float64}) = r'

#### Vector - Matrix / Matrix-Vector Products

"""
    rxp(r::Matrix{Float64}, p::Vector{Float64})

Multiply a p-vector by an r-matrix.

# Input

 - `r`     -- r-matrix
 - `p`     -- p-vector

# Output

 - `rp`    -- r * p

# Note

1) It is permissible for p and rp to be the same array.
"""
rxp(r::Matrix{Float64}, p::Vector{Float64}) = r*p

"""
    rxpv(r::Matrix{Float64}, pv::Vector{Vector{Float64}})

Multiply a pv-vector by an r-matrix.

# Input

 - `r`     -- r-matrix
 - `pv`    -- pv-vector

# Output

 - `rpv`   -- r * pv

# Note

1) The algorithm is for the simple case where the r-matrix r is not a
   function of time.  The case where r is a function of time leads to
   an additional velocity component equal to the product of the
   derivative of r and the position vector.

2) It is permissible for pv and rpv to be the same array.
"""
rxpv(r::Matrix{Float64}, pv::Vector{Vector{Float64}}) = [r*pv[1], r*pv[2]]

"""
    trxp(r::Matrix{Float64}, p::Vector{Float64})

Multiply a p-vector by the transpose of an r-matrix.

# Input

 - `r`     -- r-matrix
 - `p`     -- p-vector

# Output

 - `trp`   -- r^T * p

# Note

1) It is permissible for p and trp to be the same array.
"""
trxp(r::Matrix{Float64}, p::Vector{Float64}) = r'*p

"""
    trxpv(r::Matrix{Float64}, pv::Vector{Vector{Float64}})

Multiply a pv-vector by the transpose of an r-matrix.

# Input

 - `r`     -- r-matrix
 - `pv`    -- pv-vector

# Output

 - `trpv`  -- r^T * pv

# Note

1) The algorithm is for the simple case where the r-matrix r is not a
   function of time.  The case where r is a function of time leads to
   an additional velocity component equal to the product of the
   derivative of the transpose of r and the position vector.

2) It is permissible for pv and rpv to be the same array.
"""
trxpv(r::Matrix{Float64}, pv::Vector{Vector{Float64}}) = [r'*pv[1], r'*pv[2]]

#### Vector - Matrix / Rotation Vectors

"""
    rm2v(r::Matrix{Float64})

Express an r-matrix as an r-vector.

# Input

 - `r`     -- rotation matrix

# Output

 - `w`     -- rotation vector (Note 1)

# Note

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
    rv2m(w::Vector{Float64})

Form the r-matrix corresponding to a given r-vector.

# Input

 - `w`     -- rotation vector (Note 1)

# Output

 - `r`     -- rotation matrix

# Note

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
    pap(a::Vector{Float64}, b::Vector{Float64})

Position-angle from two p-vectors.

# Input

 - `a`     -- direction of reference point
 - `b`     -- direction of point whose PA is required

# Output

 - `θ`     -- position angle of b with respect to a (radians)

# Note

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
    pas(λa::Float64, ϕa::Float64, λb::Float64, ϕb::Float64)

Position-angle from spherical coordinates.

# Input

 - `λa`    -- longitude of point A (e.g. RA) in radians
 - `ϕa`    -- latitude of point A (e.g. Dec) in radians
 - `λb`    -- longitude of point B
 - `ϕb`    -- latitude of point B

# Output

 - `θ`     -- position angle of B with respect to A

# Note

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
    sepp(a::Vector{Float64}, b::Vector{Float64})

Angular separation between two p-vectors.

# Input

 - `a`     -- first p-vector (not necessarily unit length)
 - `b`     -- second p-vector (not necessarily unit length)

# Output

 - `θ`     -- angular separation (radians, always positive)

# Note

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
    seps(λa::Float64, ϕa::Float64, λb::Float64, ϕb::Float64)

Angular separation between two sets of spherical coordinates.

# Input

 - `λa`     -- first longitude (radians)
 - `ϕa`     -- first latitude (radians)
 - `λb`     -- second longitude (radians)
 - `ϕb`     -- second latitude (radians)

# Output

 - `θ`      -- angular separation (radians)

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
    c2s(pos::Vector{Float64})

P-vector to spherical coordinates.

# Input

 - `p`     -- p-vector

# Output

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)

# Note

1) The vector p can have any magnitude; only its direction is used.

2) If p is null, zero θ and ϕ are returned.

3) At either pole, zero θ is returned.
"""
function c2s(pos::Vector{Float64})
    NamedTuple{(:θ, :ϕ)}
    ((sum(pos[1:2].^2) == 0. ? 0. : atan(pos[2], pos[1]),
      pos[3] == 0. ? 0. : atan(pos[3], norm(pos[1:2]))))
end

"""
    p2s(pos::Vector{Float64})

P-vector to spherical polar coordinates.

# Input

 - `p`     -- p-vector

# Output

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance

# Note

1) If P is null, zero θ, ϕ and r are returned.

2) At either pole, zero θ is returned.
"""
function p2s(pos::Vector{Float64})
    @inline NamedTuple{(:θ, :ϕ, :r)}((c2s(pos)..., norm(pos)))
end

"""
    pv2s(pv::Vector{Vector{Float64}})

Convert position/velocity from Cartesian to spherical coordinates.

# Input

 - `posvel` -- position-velocity-vector

# Output

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance
 - `dθ`    -- rate of change of θ
 - `dϕ`    -- rate of change of ϕ
 - `dr`    -- rate of change of r

# Note

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
    
    x,  y,  z  = pv[norm(pv[1]) == 0.0 ? 2 : 1]
    dx, dy, dz = pv[2]

    if norm([x, y]) != 0.0
        θ  = atan(y, x)
        ϕ  = atan(z, norm([x, y]))
        dθ = (x*dy - y*dx) / (x*x+y*y)
        dϕ = (dz*(x*x+y*y) - z*(x*dx+y*dy)) / (sum([x, y, z].^2)*norm([x, y]))
    else
        θ,  ϕ  = 0.0, (z != 0.0) ? atan(z, norm([x, y])) : 0.0
        dθ, dϕ = 0.0, 0.0
    end
    r  = norm(pv[1])
    dr = norm([x, y, z]) != 0.0 ? (x*dx+y*dy+z*dz)/norm([x, y, z]) : 0.0
    (θ = θ, ϕ = ϕ, r = r, δθ = dθ, δϕ = dϕ, δr = dr)
end

"""
    s2c(θ::Float64, ϕ::Float64)

Convert spherical coordinates to Cartesian.

# Input

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)

# Output

 - `c`     -- direction cosines

"""
function s2c(θ::Float64, ϕ::Float64)
    [cos(θ)*cos(ϕ), sin(θ)*cos(ϕ), sin(ϕ)]
end

"""
    s2p(θ::Float64, ϕ::Float64, r::Float64)

Convert spherical polar coordinates to p-vector.

# Input

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance

# Output

 - `p`     -- Cartesian coordinates
"""
s2p(θ::Float64, ϕ::Float64, r::Float64) = r*s2c(θ, ϕ)

"""
    s2pv(θ::Float64, ϕ::Float64, r::Float64, dθ::Float64, dϕ::Float64,
         dr::Float64)

Convert position/velocity from spherical to Cartesian coordinates.

# Input

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance
 - `dθ`    -- rate of change of θ
 - `dϕ`    -- rate of change of ϕ
 - `dr`    -- rate of change of r

# Output

 - `pv`    -- pv-vector
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
    pdp(a::Vector{Float64}, b::Vector{Float64})

P-vector inner (=scalar=dot) product.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `r`     -- a . b
"""
pdp(a::Vector{Float64}, b::Vector{Float64}) = sum(a.*b)

"""
    pm(p::Vector{Float64}) = norm(p)

Modulus of p-vector.

# Input

 - `p`     -- p-vector

# Output

 - `r`     -- modulus
"""
pm(p::Vector{Float64}) = norm(p)

"""
    pmp(a::Vector{Float64}, b::Vector{Float64})

P-vector subtraction.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `amb`   -- a - b
"""
pmp(a::Vector{Float64}, b::Vector{Float64}) = a .- b

"""
    pn(p::Vector{Float64})

Convert a p-vector into modulus and unit vector.

# Input

 - `p`     -- p-vector

# Output

 - `r`     -- modulus
 - `u`     -- unit vector

# Note

1) If p is null, the result is null.  Otherwise the result is a unit
   vector.
"""
function pn(p::Vector{Float64})
    NamedTuple{(:modulus, :unit)}
    (norm(p) == 0 ? (0., [0., 0., 0.]) : (norm(p), p./norm(p)))
end

"""
    ppp(a::Vector{Float64}, b::Vector{Float64})

P-vector addition.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `apb`   -- a + b
"""
ppp(a::Vector{Float64}, b::Vector{Float64}) = a.+b

"""
    ppsp(a::Vector{Float64}, s::Float64, b::Vector{Float64})

P-vector plus scaled p-vector.

# Input

 - `a`     -- first p-vector
 - `s`     -- scalar (multiplier for b)
 - `b`     -- second p-vector

# Output

 - `apsb`  -- a + s*b
"""
ppsp(a::Vector{Float64}, s::Float64, b::Vector{Float64}) = a + s*b

"""
    pvdpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})

Inner (=scalar=dot) product of two pv-vectors.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `adb`   -- a . b (see note)

# Note

1) If the position and velocity components of the two pv-vectors are (
   ap, av ) and ( bp, bv ), the result, a . b, is the pair of numbers
   ( ap . bp , ap . bv + av . bp ).  The two numbers are the
   dot-product of the two p-vectors and its derivative.
"""
function pvdpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})
   [sum(a[1].*b[1]), sum(a[1].*b[2] .+ a[2].*b[1])]
end

"""
    pvm(pv::Vector{Vector{Float64}})

Modulus of pv-vector.

# Input

 - `pv`    -- pv-vector

# Output

 - `r`     -- modulus of position component
 - `s`     -- modulus of velocity component
"""
pvm(pv::Vector{Vector{Float64}}) = sqrt.(sum.([pv[1].^2, pv[2].^2]))

"""
    pvmpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})

Subtract one pv-vector from another.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `amb`   -- a - b
"""
pvmpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}}) = [a[1] .- b[1], a[2] .- b[2]]

"""
    pvppv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})

Add one pv-vector to another.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `apb`   -- a + b
"""
pvppv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}}) = [a[1] .+ b[1], a[2] .+ b[2]]

"""
    pvu(dt::Float64, pv::Vector{Vector{Float64}})

Update a pv-vector.

# Input

 - `dt`    -- time interval
 - `pv`    -- pv-vector

# Output

 - `upv`   -- p updated, v unchanged

# Note

1) "Update" means "refer the position component of the vector to a new
   date dt time units from the existing date".

2) The time units of dt must match those of the velocity.
"""
pvu(dt::Float64, pv::Vector{Vector{Float64}}) = [pv[1] .+ dt.*pv[2], pv[2]]

"""
    pvup(dt::Float64, pv::Vector{Vector{Float64}})

Update a pv-vector, discarding the velocity component.

# Input

 - `dt`    -- time interval
 - `pv`    -- pv-vector

# Output

 - `p`     -- p-vector

# Note

1) "Update" means "refer the position component of the vector to a new
   date dt time units from the existing date".

2) The time units of dt must match those of the velocity.
"""
pvup(dt::Float64, pv::Vector{Vector{Float64}}) = pv[1] .+ dt*pv[2]

"""
    pvxpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})

Outer (=vector=cross) product of two pv-vectors.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `axb`   -- a x b

# Note

1) If the position and velocity components of the two pv-vectors are (
   ap, av ) and ( bp, bv ), the result, a x b, is the pair of vectors
   ( ap x bp, ap x bv + av x bp ).  The two vectors are the
   cross-product of the two p-vectors and its derivative.
"""
function pvxpv(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})
    [vec2mat(a[1])*b[1], vec2mat(a[1])*b[2] .+ vec2mat(a[2])*b[1]]
end

"""
    pxp(a::Vector{Float64}, b::Vector{Float64})

P-vector outer (=vector=cross) product.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `axb`   -- a x b
"""
pxp(a::Vector{Float64}, b::Vector{Float64}) = vec2mat(a)*b

"""
    s2xpv(s1::Float64, s2::Float64, pv::Vector{Vector{Float64}})

Multiply a pv-vector by two scalars.

# Input

 - `s1`    -- scalar to multiply position component by
 - `s2`    -- scalar to multiply velocity component by
 - `pv`    -- pv-vector

# Output

 - `spv`   -- pv-vector: p scaled by s1, v scaled by s2
"""
s2xpv(s1::Float64, s2::Float64, pv::Vector{Vector{Float64}}) = [s1*pv[1], s2*pv[2]]

"""
    sxp(s::Float64, p::Vector{Float64})

Multiply a p-vector by a scalar.

# Input

 - `s`     -- scalar
 - `p`     -- p-vector

# Output

 - `sp`    -- s * p
"""
sxp(s::Float64, p::Vector{Float64}) = s*p

"""
    sxpv(s::Float64, pv::Vector{Vector{Float64}})

Multiply a pv-vector by a scalar.

# Input

 - `s`     -- scalar
 - `pv`    -- pv-vector

# Output

 - `spv`   -- s * pv
"""
sxpv(s::Float64, pv::Vector{Vector{Float64}}) = [s*pv[1], s*pv[2]]
