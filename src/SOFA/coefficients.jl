####    Astronomy / Fundamental Arguments    ####

""" 
    fad03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean elongation of the
Moon from the Sun.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `D`     -- D, radians (Note 2)

# Note

1) Though Δt is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and is
   from Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
"""
fad03(Δt::Float64) = deg2rad(1/3600)*rem(Polynomial(D_2003A, :Δt)(Δt), ARCSECPER2PI)

"""
    fae03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of
Earth.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Earth, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111
"""
fae03(Δt::Float64) = mod2pi(Polynomial(lea_2003, :Δt)(Δt))

"""
    faf03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of the
Moon minus mean longitude of the ascending node.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `F`     -- F, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and is
   from Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
"""
faf03(Δt::Float64) = deg2rad(1/3600)*rem(Polynomial(F_2003A, :Δt)(Δt), ARCSECPER2PI)

"""
    faju03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of
Jupiter.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Jupiter, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111
"""
faju03(Δt::Float64) = mod2pi(Polynomial(lju_2003, :Δt)(Δt))

"""
    fal03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean anomaly of the
Moon.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `l`     -- l, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and is
   from Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
"""
fal03(Δt::Float64) = deg2rad(rem(Polynomial(l0_2003A, :Δt)(Δt), ARCSECPER2PI)/3600)

"""
    falp03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean anomaly of the
Sun.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output\

 - `l`     -- l', radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and is
   from Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
"""
falp03(Δt::Float64) = deg2rad(1/3600)*rem(Polynomial(l1_2003A, :Δt)(Δt), ARCSECPER2PI)

"""
    fama03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of Mars.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Mars, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111
"""
fama03(Δt::Float64) = mod2pi(Polynomial(lma_2003, :Δt)(Δt))

"""
    fame03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of
Mercury.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Mercury, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111
"""
fame03(Δt::Float64) = mod2pi(Polynomial(lme_2003, :Δt)(Δt))

"""
    fane03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of
Neptune.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Neptune, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and is
   adapted from Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
"""
fane03(Δt::Float64) = mod2pi(Polynomial(lne_2003, :Δt)(Δt))

"""
    faom03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of the
Moon's ascending node.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output
 - `Ω`     -- Ω, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and is
   from Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J., 1994, Astron.Astrophys. 282, 663-683.
"""
faom03(Δt::Float64) = deg2rad(1/3600)*rem(Polynomial(Ω_2003A, :Δt)(Δt), ARCSECPER2PI)

"""
    fapa03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): general accumulated
precession in longitude.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `long`  -- general precession in longitude, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003).  It
   is taken from Kinoshita & Souchay (1990) and comes originally from
   Lieske et al. (1977).

# References

Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.  48,
187

Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
Astron.Astrophys. 58, 1-16

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
fapa03(Δt::Float64) = mod2pi(Polynomial(lge_2003, :Δt)(Δt))

"""
    fasa03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of
Saturn.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Saturn, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111
"""
fasa03(Δt::Float64) = mod2pi(Polynomial(lsa_2003, :Δt)(Δt))

"""
    faur03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of
Uranus.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Uranus, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and is
   adapted from Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
"""
faur03(Δt::Float64) = mod2pi(Polynomial(lur_2003, :Δt)(Δt))

"""
    fave03(Δt::Float64)

Fundamental argument, IERS Conventions (2003): mean longitude of
Venus.

# Input

 - `Δt`    -- TDB, Julian centuries since J2000.0 (Note 1)

# Output

 - `mlon`  -- mean longitude of Venus, radians (Note 2)

# Note

1) Though t is strictly TDB, it is usually more convenient to use TT,
   which makes no significant difference.

2) The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111
"""
fave03(Δt::Float64) = mod2pi(Polynomial(lve_2003, :Δt)(Δt))
