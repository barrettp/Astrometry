module Astrometry

export JD2000, MJD0, calendar_mjd
export iau_1980_nutation, iau_1980_obliquity, iau_1976_precession
export iau_2000_equinox_complement, iau_2000a_nutation, iau_2000b_nutation
export iau_2000_position, iau_2000_precession
export iau_2006a_nutation, iau_2006_cio_locator, iau_2006_cip_xy, iau_2006_obliquity
export iau_2006_precession, iau_2006_tdb_tt
export equinox, precession_nutation, proper_motion, radec2rad
export SOFA

using Polynomials

include("constants.jl")
include("constants2011.jl")
include("constantsplanet.jl")
include("util.jl")
include("model1980.jl")
include("model2000.jl")
include("model2006.jl")
include("astrom.jl")
include("earth.jl")

"""
Standards of Fundamental Astronomy (SOFA)

Version: 2021-05-21
"""
module SOFA

#   Astronomy: Calenders
export cal2jd, epb, epb2jd, epj, epj2jd, jd2cal, jdcalf
#   Astronomy: Astrometry
export ab, apcg, apcg13, apci, apci13, apco, apco13, apcs, apcs13, aper, aper13,
    apio, apio13, atcc13, atccq, atci13, atciq, atciqn, atciqz, atco13, atic13,
    aticq, aticqn, atio13, atioq, atoc13, atoi13, atoiq, ld, ldn, ldsun, pmpx,
    pmsafe, pvtob, refco
#   Astronomy: Ephemerides
export epv00, moon98, plan94
#   Astronomy: Fundamental Constants
export fad03, fae03, faf03, faju03, fal03, falp03, fama03, fame03, fane03,
    faom03, fapa03, fasa03, faur03, fave03
#   Astronomy: Precession-Nutation-Polar Motion
export bi00, bp00, bp06, bpn2xy, c2i00a, c2i00b, c2i06a, c2ibpn, c2ixy, c2ixys,
    c2t00a, c2t00b, c2t06a, c2tcio, c2teqx, c2tpe, c2txy, eo06a, eors, fw2m, f2wxy,
    ltp, ltpb, ltpecl, ltpequ, num00a, num00b, num06a, nut00a, nut00b, nut06a, nut80,
    nutm80, obl06, obl80, po6e, pb06, pfw06, pmat00, pmat06, pmat76, pn00, pn00a,
    pn00b, pn06, pn06a, pnm00a, pnm00b, pnm06a, pnm80, pom00, pr00, prec76, s00,
    s00a, s00b, s06, s06a, sp00, xy06, xys00a, xys00b, xys06a
#   Astronomy: Rotation and Time
export ee00, ee00a, ee00b, ee06a, eect00, eqeq94, era00, gmst00, gmst06, gmst82,
    gst00a, gst00b, gst06, gst06a, gst94
#   Astronomy: Space Motion
export pvstar, starpv
#   Astronomy: Star Catalogs
export fk425, fk45z, fk524, fk52h, fk54z, fk5hip, fk5hz, h2fk5, hfk5z, starpm
#   Astronomy: Ecliptic Coordinates
export eceq06, ecm06, eqec06, lteceq, ltecm, lteqec
#   Astronomy: Galactic Coordinates
export g2icrs, icrs2g
#   Astronomy: Geodetic-Geocentric Coordinates
export eform, gc2gd, gc2gde, gd2gc, gd2gce
#   Astronomy: Timescales
export d2dtf, dat, dtdb, dtf2d, taitt, taiut1, taiutc, tcbtdb, tcgtt, tdbtcb,
    tdbtt, tttai, tttcg, tttdb, ttut1, ut1tai, ut1tt, ut1utc, utctai, utcut1
#   Astronomy: Horizontal-Equatorial
export ae2hd, hd2ae, hd2pa
#   Astronomy: Gnomonic
export tpors, tporv, tpsts, tpstv, tpxes, tpxev
#   Vector-Matrix: Angle Operations
export a2af, a2tf, af2a, anp, anpm, d2tf, tf2a, tf2d
#   Vector-Matrix: Build Rotations
export rx, ry, rz
#   Vector-Matrix: Copy-Extend-Extract
export cp, cpv, cr, p2pv, pv2p
#   Vector-Matrix: Initialization
export ir, zp, zpv, zr
#   Vector-Matrix: Matrix Operations
export rxr, tr
#   Vector-Matrix: Matrix-Vector Products
export rxp, rxpv, trxp, trxpv
#   Vector-Matrix: Roation Vectors
export rm2v, rv2m
#   Vector-Matrix: Separation Vectors and Angles
export pap, pas, sepp, seps
#   Vector-Matrix: Spherical-Cartesian
export c2s, p2s, pv2s, s2c, s2p, s2pv
#   Vector-Matrix: Vector Operations
export pdp, pm, pmp, pn, ppp, ppsp, pvdpv, pvm, pvmpv, pvppv, pvu, pvup, pvxpv,
    pxp, s2xpv, sxp, sxpv

import ..ARCSECPER2PI, ..ASECPERRAD, ..ASTRUNIT, ..AULIGHT, ..LIGHTSPEED, ..SCHWARZRADIUS
import ..SECPERDAY, ..DAYPERYEAR, ..DAYINYEAR1900, ..DAYINYEAR2000, ..MODJULDAY0, ..JDMIN, ..JDMAX
import ..ELB, ..ELG, ..TT_MINUS_TAI, ..TDB0, ..BD1900, ..MJD77, ..MJD0, ..MJD00, ..MJDAY0, ..JD2000
import ..JULIANDAY2000, ..GK, ..COSEPS, ..SINEPS
import ..wgs72_radius, ..wgs72_oblate, ..grs80_radius, ..grs80_oblate, ..wgs84_radius,
    ..wgs84_oblate
#  1976 precession model
import ..ζT_1976, ..θT_1976, ..ζA_1976, ..zA_1976, ..θA_1976
import ..ψ_1977, ..ω_1977, ..χ_1977
#  1980 nutation model
import ..l0_1980t, ..l1_1980t, ..F_1980t, ..D_1980t, ..Ω_1980t, ..l0_1980, ..l1_1980,
    ..F_1980, ..D_1980, ..Ω_1980, ..ϵ_1980, ..iau_1980_nutation_series
#  1982 Greenwich mean sidereal time
import ..gmst_1982
#  1994 model
import ..l_1994, ..equinox_1994
#
import ..ϵsun_1994, ..ϵmsun_1994, ..D_1994, ..ϵju_1994, ..ϵsa_1994, ..topo_1994, ..λmoon_1994
import ..mass_1994, ..a_1994, ..λ_1994, ..e_1994, ..π_1994, ..i_1994, ..ω_1994,
    ..p_1994, ..q_1994, ..a_cos_1994, ..a_sin_1994, ..λ_cos_1994, ..λ_sin_1994
import ..mass_plan_1994_0, ..mass_plan_1994_2
import ..dmoon_1998, ..lsun_1998, ..lmoon_1998, ..fmoon_1998, ..lr_1998, ..b_1998
import ..a_1, ..a_2, ..a_3, ..a_l, ..a_b, ..r0, ..efac
#  2000 model
import ..era_2000, ..gmst_2000, ..tio_2000, ..ϵ0_2000, ..p03_2000
import ..iau_2000_equinox_0_series, ..iau_2000_equinox_1_series
import ..ψ_bias_2000, ..ϵ_bias_2000, ..ψ_corr_2000, ..ϵ_corr_2000, ..icrs_ra_2000,
    ..j2_corr_2000
import ..r_gal_icrs, ..iau_2000_bcrs
import ..l1_2000A, ..D_2000A, ..l0_2000A_planet, ..F_2000A_planet, ..D_2000A_planet,
    ..Ω_2000A_planet
#  2000A nutation model
import ..iau_2000A_nutation_lunisolar_series, ..iau_2000A_nutation_planetary_series
import ..sp_2000A, ..s0_2000A, ..s1_2000A, ..s2_2000A, ..s3_2000A, ..s4_2000A
#  2000B nutation model
import ..l1_2000B, ..D_2000B, ..l0_2000B, ..F_2000B, ..Ω_2000B, ..ψ_2000B_planet, ..ϵ_2000B_planet
import ..iau_2000B_nutation_lunisolar_series
import ..sun_earth_x_0, ..sun_earth_x_1, ..sun_earth_x_2, ..sun_earth_y_0, ..sun_earth_y_1,
    ..sun_earth_y_2, ..sun_earth_z_0, ..sun_earth_z_1, ..sun_earth_z_2, ..sun_earth_x_0,
    ..sun_earth_x_1, ..sun_earth_x_2, ..sun_earth_y_0, ..sun_earth_y_1, ..sun_earth_y_2,
    ..sun_earth_z_0, ..sun_earth_z_1, ..sun_earth_z_2
import ..bary_sun_x_0, ..bary_sun_x_1, ..bary_sun_x_2, ..bary_sun_y_0, ..bary_sun_y_1,
    ..bary_sun_y_2, ..bary_sun_z_0, ..bary_sun_z_1, ..bary_sun_z_2, ..bary_sun_x_0, ..bary_sun_x_1,
    ..bary_sun_x_2, ..bary_sun_y_0, ..bary_sun_y_1, ..bary_sun_y_2, ..bary_sun_z_0, ..bary_sun_z_1,
    ..bary_sun_z_2
import ..D_2003A, ..lea_2003, ..F_2003A, ..lju_2003, ..l0_2003A, ..l1_2003A, ..lma_2003,
    ..lme_2003, ..lne_2003, ..Ω_2003A, ..lge_2003, ..lsa_2003, ..lur_2003, ..lve_2003,
    ..lne_2003mhb
import ..tdb_tt_2003_0, ..tdb_tt_2003_1, ..tdb_tt_2003_2, ..tdb_tt_2003_3, ..tdb_tt_2003_4
#  2000A precession-nutation model
import ..ϵ0_2006, ..ψA_2006, ..ωA_2006, ..PA_2006, ..QA_2006, ..πA_2006, ..ΠA_2006, ..ϵA_2006,
    ..χA_2006, ..ζA_2006, ..θA_2006, ..zA_2006, ..pA_2006, ..γF_2006, ..ϕF_2006, ..ψF_2006
#  2000B precession-nutation model
import ..ϵB_2006, ..γB_2006, ..ϕB_2006, ..ψB_2006
import ..cio_s_2006, ..iau_2006_equinox_0_series, ..iau_2006_equinox_1_series,
    ..iau_2006_equinox_2_series, ..iau_2006_equinox_3_series, ..iau_2006_equinox_4_series
#  2006 Celestial Intermediate Pole
import ..cip_x_2006, ..cip_y_2006, ..cip_lunisolar_2006, ..cip_planetary_2006,
    ..cip_pointer_2006, ..cip_amplitude_2006
#  2006 Earth rotation model
import ..Ω_Earth_2003
import ..gmst_2006
import ..ϵ0_2010, ..η0_2010, ..dα0_2010
import ..ecl_ϕ_2011, ..ecl_pA_0_2011, ..ecl_qA_0_2011, ..ecl_pA_c_2011, ..ecl_pA_s_2011,
    ..ecl_qA_c_2011, ..ecl_qA_s_2011
import ..equ_ϕ_2011, ..equ_xA_0_2011, ..equ_yA_0_2011, ..equ_xA_c_2011, ..equ_xA_s_2011,
    ..equ_yA_c_2011, ..equ_yA_s_2011
import ..calendar2MJD, ..Rx, ..Ry, ..Rz, ..vec2mat, ..ephem_position, ..ephem_velocity

using LinearAlgebra, Polynomials

include("sofa.jl")

end

end

#   Issues
#   1. Timing accuracy, i.e. significant digits.
