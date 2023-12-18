#### Astronomy / Gnomonic

"""
    tpors(ξ::Float64, η::Float64, a::Float64, b::Float64)

In the tangent plane projection, given the rectangular coordinates of
a star and its spherical coordinates, determine the spherical
coordinates of the tangent point.

# Input

 - `ξ, η`  -- rectangular coordinates of star image (Note 2)
 - `a, b`  -- star's spherical coordinates (Note 3)

# Output

 - `a01, b01` -- tangent point's spherical coordinates, Soln. 1
 - `a02, b02` -- tangent point's spherical coordinates, Soln. 2

# Note

1) The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2) The eta axis points due north in the adopted coordinate system.  If
   the spherical coordinates are observed (RA,Dec), the tangent plane
   coordinates (ξ,η) are conventionally called the "standard
   coordinates".  If the spherical coordinates are with respect to a
   right-handed triad, (ξ,η) are also right-handed.  The units of
   (ξ,η) are, effectively, radians at the tangent point.

3) All angular arguments are in radians.

4) The angles a01 and a02 are returned in the range 0-2π.  The angles
   b01 and b02 are returned in the range +/-π, but in the usual,
   non-pole-crossing, case, the range is +/-π/2.

5) Cases where there is no solution can arise only near the poles.
   For example, it is clearly impossible for a star at the pole itself
   to have a non-zero ξ value, and hence it is meaningless to ask
   where the tangent point would have to be to bring about this
   combination of ξ and dec.

6) Also near the poles, cases can arise where there are two useful
   solutions.  The return value indicates whether the second of the
   two solutions returned is useful; 1 indicates only one useful
   solution, the usual case.

7) The basis of the algorithm is to solve the spherical triangle PSC,
   where P is the north celestial pole, S is the star and C is the
   tangent point.  The spherical coordinates of the tangent point are
   [a0,b0]; writing ρ^2 = (ξ^2+η^2) and r^2 = (1+ρ^2), side c is then
   (π/2-b), side p is sqrt(ξ^2+η^2) and side s (to be found) is
   (π/2-b0).  Angle C is given by sin(C) = ξ/ρ and cos(C) = η/ρ.
   Angle P (to be found) is the longitude difference between star and
   tangent point (a-a0).

8) This function is a member of the following set:

       spherical      vector         solve for

       eraTpxes      eraTpxev          ξ,η
       eraTpsts      eraTpstv          star
     > eraTpors <    eraTporv         origin

# References

Calabretta M.R. & Greisen, E.W., 2002, "Representations of celestial
coordinates in FITS", Astron.Astrophys. 395, 1077

Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987,
Chapter 13.
"""
function tpors(ξ::Float64, η::Float64, a::Float64, b::Float64)
    r  = sqrt(1.0 + η*η + ξ*ξ)
    w2 = (r*cos(b))^2 - ξ*ξ
    if w2 >= 0.0
        wp = ξ == 0.0 && sqrt(w2) == 0.0 ? 1.0 : sqrt(w2)
        a01, b01 = mod2pi(a - atan(ξ, wp)), atan(r*sin(b) - η*wp, r*sin(b)*η + wp)
        wn = -wp
        a02, b02 = mod2pi(a - atan(ξ, wn)), atan(r*sin(b) - η*wn, r*sin(b)*η + wn)
        res = abs(r*sin(b)) < 1.0 ? (a01, b01, nothing, nothing) : (a01, b01, a02, b02)
    else
        res = (nothing, nothing, nothing, nothing)
    end
    NamedTuple{(:a01, :b01, :a02, :b02)}(res)
end

"""
    tporv(ξ::Float64, η::Float64, v::Vector{Float64})

In the tangent plane projection, given the rectangular coordinates of
a star and its direction cosines, determine the direction cosines of
the tangent point.

# Input

 - `ξ, η`  -- rectangular coordinates of star image (Note 2)
 - `v`     -- star's direction cosines (Note 3)

# Output

 - `v01`   -- tangent point's direction cosines, Solution 1
 - `v02`   -- tangent point's direction cosines, Solution 2

# Note

1) The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2) The eta axis points due north in the adopted coordinate system.  If
   the direction cosines represent observed (RA,Dec), the tangent
   plane coordinates (ξ,η) are conventionally called the "standard
   coordinates".  If the direction cosines are with respect to a
   right-handed triad, (ξ,η) are also right-handed.  The units of
   (ξ,η) are, effectively, radians at the tangent point.

3) The vector v must be of unit length or the result will be wrong.

4) Cases where there is no solution can arise only near the poles.
   For example, it is clearly impossible for a star at the pole itself
   to have a non-zero xi value, and hence it is meaningless to ask
   where the tangent point would have to be.

5) Also near the poles, cases can arise where there are two useful
   solutions.  The return value indicates whether the second of the
   two solutions returned is useful; 1 indicates only one useful
   solution, the usual case.

6) The basis of the algorithm is to solve the spherical triangle PSC,
   where P is the north celestial pole, S is the star and C is the
   tangent point.  Calling the celestial spherical coordinates of the
   star and tangent point (a,b) and (a0,b0) respectively, and writing
   ρ^2 = (ξ^2+η^2) and r^2 = (1+ρ^2), and transforming the vector v
   into (a,b) in the normal way, side c is then (π/2-b), side p is
   sqrt(ξ^2+η^2) and side s (to be found) is (π/2-b0), while angle C
   is given by sin(C) = ξ/ρ and cos(C) = η/ρ; angle P (to be found) is
   (a-a0).  After solving the spherical triangle, the result (a0,b0)
   can be expressed in vector form as v0.

7) This function is a member of the following set:

       spherical      vector         solve for

       eraTpxes      eraTpxev         xi,eta
       eraTpsts      eraTpstv          star
       eraTpors    > eraTporv <       origin

# References

Calabretta M.R. & Greisen, E.W., 2002, "Representations of celestial
coordinates in FITS", Astron.Astrophys. 395, 1077

Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987,
Chapter 13.
"""
function tporv(ξ::Float64, η::Float64, v::Vector{Float64})
    r  = sqrt(1.0 + η*η + ξ*ξ)
    w2 = r*r*sum(v[1:2].^2) - ξ*ξ
    if w2 > 0.0
        wp = sqrt(w2)
        cp = (r*v[3]*η + wp)/((1.0+η*η)*sqrt(r*r*sum(v[1:2].^2)^2))
        v01 = [cp*(v[1]*wp+v[2]*ξ), cp*(v[2]*wp-v[1]*ξ), (r*v[3]-η*wp)/(1.0+η*η)]
        wn = -sqrt(w2)
        cn = (r*v[3]*η + wn)/((1.0+η*η)*sqrt(r*r*sum(v[1:2].^2)^2))
        v02 = [cn*(v[1]*wn+v[2]*ξ), cn*(v[2]*wn-v[1]*ξ), (r*v[3]-η*wn)/(1.0+η*η)]
        res = abs(r*v[3]) < 1.0 ? (v01, nothing) : (v01, v02)
    else
        res = (nothing, nothing)
    end
    NamedTuple{(:v01, :v02)}(res)
end

"""
    tpsts(ξ::Float64, η::Float64, a0::Float64, b0::Float64)

In the tangent plane projection, given the star's rectangular
coordinates and the spherical coordinates of the tangent point, solve
for the spherical coordinates of the star.

# Input

 - `ξ, η`   -- rectangular coordinates of star image (Note 2)
 - `a0, b0` -- tangent point's spherical coordinates

# Output

 - `a, b`   -- star's spherical coordinates

# Note

1) The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2) The eta axis points due north in the adopted coordinate system.  If
   the spherical coordinates are observed (RA,Dec), the tangent plane
   coordinates (ξ,η) are conventionally called the "standard
   coordinates".  If the spherical coordinates are with respect to a
   right-handed triad, (ξ,η) are also right-handed.  The units of
   (ξ,η) are, effectively, radians at the tangent point.

3) All angular arguments are in radians.

4) This function is a member of the following set:

       spherical      vector         solve for

       eraTpxes      eraTpxev         xi,eta
     > eraTpsts <    eraTpstv          star
       eraTpors      eraTporv         origin

# References

Calabretta M.R. & Greisen, E.W., 2002, "Representations of celestial
coordinates in FITS", Astron.Astrophys. 395, 1077

Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987,
Chapter 13.
"""
function tpsts(ξ::Float64, η::Float64, a0::Float64, b0::Float64)
    NamedTuple{(:a, :b)}
    ((mod2pi(atan(ξ, cos(b0) - η*sin(b0)) + a0),
      atan(sin(b0) + η*cos(b0), sqrt(ξ*ξ+(cos(b0) - η*sin(b0))^2))))
end

"""
    tpstv(ξ::Float64, η::Float64, v0::Vector{Float64})

In the tangent plane projection, given the star's rectangular
coordinates and the direction cosines of the tangent point, solve for
the direction cosines of the star.

# Input

 - `ξ, η`  -- rectangular coordinates of star image (Note 2)
 - `v0`    -- tangent point's direction cosines

# Output

 - `v`     -- star's direction cosines

# Note

1) The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2) The η axis points due north in the adopted coordinate system.  If
   the direction cosines represent observed (RA,Dec), the tangent
   plane coordinates (ξ,η) are conventionally called the "standard
   coordinates".  If the direction cosines are with respect to a
   right-handed triad, (ξ,η) are also right-handed.  The units of
   (ξ,η) are, effectively, radians at the tangent point.

3) The method used is to complete the star vector in the (xi,eta)
   based triad and normalize it, then rotate the triad to put the
   tangent point at the pole with the x-axis aligned to zero
   longitude.  Writing (a0,b0) for the celestial spherical coordinates
   of the tangent point, the sequence of rotations is (b-π/2) around
   the x-axis followed by (-a-π/2) around the z-axis.

4) If vector v0 is not of unit length, the returned vector v will
   be wrong.

5) If vector v0 points at a pole, the returned vector v will be
   based on the arbitrary assumption that the longitude coordinate
   of the tangent point is zero.

6) This function is a member of the following set:

       spherical      vector         solve for

       eraTpxes      eraTpxev           ξ,η
       eraTpsts    > eraTpstv <        star
       eraTpors      eraTporv         origin

# References

Calabretta M.R. & Greisen, E.W., 2002, "Representations of celestial
coordinates in FITS", Astron.Astrophys. 395, 1077

Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987,
Chapter 13.
"""
function tpstv(ξ::Float64, η::Float64, v0::Vector{Float64})
    x, y, z = v0[1:3]
    if sqrt(sum(v0[1:2].^2)) == 0.0
        r = x = 1e-20
    else
        r = sqrt(sum(v0[1:2].^2))
    end
    f = sqrt(1.0 + ξ*ξ + η*η)
    [x - (ξ*y + η*x*z)/r, y + (ξ*x - η*y*z)/r, z + η*r]/f
end

"""
    tpxes(a::Float64, b::Float64, a0::Float64, b0::Float64)

In the tangent plane projection, given celestial spherical coordinates
for a star and the tangent point, solve for the star's rectangular
coordinates in the tangent plane.

# Input

 - `a, b`   -- star's spherical coordinates
 - `a0, b0` -- tangent point's spherical coordinates

# Output

 - `ξ, η`   -- rectangular coordinates of star image (Note 2)

# Note

1) The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2) The η axis points due north in the adopted coordinate system.  If
   the spherical coordinates are observed (RA,Dec), the tangent plane
   coordinates (ξ,η) are conventionally called the "standard
   coordinates".  For right-handed spherical coordinates, (ξ,η) are
   also right-handed.  The units of (ξ,η) are, effectively, radians at
   the tangent point.

3) All angular arguments are in radians.

4) This function is a member of the following set:

       spherical      vector         solve for

     > eraTpxes <    eraTpxev           ξ,η
       eraTpsts      eraTpstv          star
       eraTpors      eraTporv         origin

# References

Calabretta M.R. & Greisen, E.W., 2002, "Representations of celestial
coordinates in FITS", Astron.Astrophys. 395, 1077

Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987,
Chapter 13.
"""
function tpxes(a::Float64, b::Float64, a0::Float64, b0::Float64)
    d = (sin(b)*sin(b0) + cos(b)*cos(b0)*cos(a - a0))
    if d > TINY
        nothing
    elseif d >= 0.0
        d = TINY
    elseif d > -TINY
        d = -TINY
    else
        nothing
    end
    NamedTuple{(:ξ, :η)}
    (cos(b)*sin(a - a0)/d, (sin(b)*cos(b0) - cos(b)*sin(b0)*cos(a - a0))/d)
end

"""
    tpxev(v::Vector{Float64}, v0::Vector{Float64})

In the tangent plane projection, given celestial direction cosines for
a star and the tangent point, solve for the star's rectangular
coordinates in the tangent plane.

# Input

 - `v`     -- direction cosines of star (Note 4)
 - `v0`    -- direction cosines of tangent point (Note 4)

# Output

 - `ξ, η`  -- tangent plane coordinates of star

# Note

1) The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2) The eta axis points due north in the adopted coordinate system.  If
   the direction cosines represent observed (RA,Dec), the tangent
   plane coordinates (ξ,η) are conventionally called the "standard
   coordinates".  If the direction cosines are with respect to a
   right-handed triad, (ξ,η) are also right-handed.  The units of
   (ξ,η) are, effectively, radians at the tangent point.

3) The method used is to extend the star vector to the tangent plane
   and then rotate the triad so that (x,y) becomes (ξ,η).  Writing
   (a,b) for the celestial spherical coordinates of the star, the
   sequence of rotations is (a+π/2) around the z-axis followed by
   (π/2-b) around the x-axis.

4) If vector v0 is not of unit length, or if vector v is of zero
   length, the results will be wrong.

5) If v0 points at a pole, the returned (ξ,η) will be based on the
   arbitrary assumption that the longitude coordinate of the tangent
   point is zero.

6) This function is a member of the following set:

       spherical      vector         solve for

       eraTpxes    > eraTpxev <         ξ,η
       eraTpsts      eraTpstv          star
       eraTpors      eraTporv         origin

# References

Calabretta M.R. & Greisen, E.W., 2002, "Representations of celestial
coordinates in FITS", Astron.Astrophys. 395, 1077

Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987,
Chapter 13.
"""
function tpxev(v::Vector{Float64}, v0::Vector{Float64})
    x,   y,  z = v[1:3]
    x0, y0, z0 = v0[1:3]
    if sqrt(sum(v0[1:2].^2)) == 0.0
        r = x0 = 1e-20
    else
        r = sqrt(sum(v0[1:2].^2))
    end
    d = x*x0 + y*y0 + z*z0
    if d > TINY
        nothing
    elseif d >= 0.0
        d = TINY
    elseif d > -TINY
        d = -TINY
    else
        nothing
    end
    NamedTuple{(:ξ, :η)}((y*x0 - x*y0, z*sum(v0[1:2].^2) - z0*(x*x0 + y*y0))./(d*r))
end
