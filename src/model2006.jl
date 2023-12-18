#   IAU 2006 Model
#=
struct S06coef
    fa::Vector{Int8}
    sn::Float64
    cs::Float64
end
=#

include("constants2006.jl")
include("constants2006F.jl")

"""
    iau_2006_cio_locator(date, coords)

Locate the Celestial Intermediate Origin locator, s, on the equator of the
Celestial Intermediate Pole, given the CIP X and Y coordinates. (IAU 2000A and 2006).
"""
function iau_2006_cio_locator(date, coord)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ϕ = deg2rad.(rem.([Polynomial(l0_2003A, :Δt)(Δt),
        Polynomial(l1_2003A, :Δt)(Δt), Polynomial( F_2003A, :Δt)(Δt),
        Polynomial( D_2003A, :Δt)(Δt), Polynomial( Ω_2003A, :Δt)(Δt)], ARCSECPER2PI)/3600)
    append!(ϕ, [mod2pi(Polynomial(lve_2003, :Δt)(Δt)),
                mod2pi(Polynomial(lea_2003, :Δt)(Δt)), Polynomial(lge_2003, :Δt)(Δt)])

    ϕ0 = vcat([t.n' for t in iau_2006_equinox_0_series]...)*ϕ
    a0 = vcat([t.a' for t in iau_2006_equinox_0_series]...)
    ϕ1 = vcat([t.n' for t in iau_2006_equinox_1_series]...)*ϕ
    a1 = vcat([t.a' for t in iau_2006_equinox_1_series]...)
    ϕ2 = vcat([t.n' for t in iau_2006_equinox_2_series]...)*ϕ
    a2 = vcat([t.a' for t in iau_2006_equinox_2_series]...)
    ϕ3 = vcat([t.n' for t in iau_2006_equinox_3_series]...)*ϕ
    a3 = vcat([t.a' for t in iau_2006_equinox_3_series]...)
    ϕ4 = vcat([t.n' for t in iau_2006_equinox_4_series]...)*ϕ
    a4 = vcat([t.a' for t in iau_2006_equinox_4_series]...)

    deg2rad(Polynomial(cio_s_2006 .+ [
        sum(a0[:,1].*sin.(ϕ0) .+ a0[:,2].*cos.(ϕ0)),
        sum(a1[:,1].*sin.(ϕ1) .+ a1[:,2].*cos.(ϕ1)),
        sum(a2[:,1].*sin.(ϕ2) .+ a2[:,2].*cos.(ϕ2)),
        sum(a3[:,1].*sin.(ϕ3) .+ a3[:,2].*cos.(ϕ3)),
        sum(a4[:,1].*sin.(ϕ4) .+ a4[:,2].*cos.(ϕ4)), 0.0], :Δt)(Δt)/3600) -
            coord[1]*coord[2]/2
end

function iau_2006_cip_xy(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    #  lunar, solar, and planetary longitudes
    ϕ = deg2rad(1/3600).*rem.([
        Polynomial(l0_2003A, :Δt)(Δt),
        Polynomial(l1_2003A, :Δt)(Δt),
        Polynomial( F_2003A, :Δt)(Δt),
        Polynomial( D_2003A, :Δt)(Δt),
        Polynomial( Ω_2003A, :Δt)(Δt)], ARCSECPER2PI)
    append!(ϕ, [
        Polynomial( lme_2003, :Δt)(Δt),
        Polynomial( lve_2003, :Δt)(Δt),
        Polynomial( lea_2003, :Δt)(Δt),
        Polynomial( lma_2003, :Δt)(Δt),
        Polynomial( lju_2003, :Δt)(Δt),
        Polynomial( lsa_2003, :Δt)(Δt),
        Polynomial( lur_2003, :Δt)(Δt),
        Polynomial( lne_2003, :Δt)(Δt),
        Polynomial( lge_2003, :Δt)(Δt)])
    
    #  polynomial part of precession-nutation
    xypr = [sum(Polynomial(cip_x_2006, :Δt)(Δt)),
            sum(Polynomial(cip_y_2006, :Δt)(Δt))]

    # !!! The following code can be improved by rearranging the data arrays.

    jaxy = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1] .+ 1
    jasc = [0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0] .+ 1
    japt = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]

    #  nutation periodic terms, planetary
    xypl = [0.0, 0.0]
    ialast = length(cip_amplitude_2006)
    for ifreq = length(cip_planetary_2006):-1:1
        sc = sincos(sum(cip_planetary_2006[ifreq].*ϕ))
        ia = cip_pointer_2006[ifreq + length(cip_lunisolar_2006)]
        for i = ialast:-1:ia
            xypl[jaxy[i-ia+1]] += cip_amplitude_2006[i] * sc[jasc[i-ia+1]] * Δt^japt[i-ia+1]
        end
        ialast = ia-1
    end

    #  nutation periodic terms, luni-solar
    xyls = [0.0, 0.0]
    for ifreq = length(cip_lunisolar_2006):-1:1
        sc = sincos(sum(cip_lunisolar_2006[ifreq].*ϕ[1:5]))
        ia = cip_pointer_2006[ifreq]
        for i = ialast:-1:ia
            xyls[jaxy[i-ia+1]] += cip_amplitude_2006[i] * sc[jasc[i-ia+1]] * Δt^japt[i-ia+1]
        end
        ialast = ia-1
    end

    deg2rad.((xypr .+ (xyls .+ xypl)./1e6)/3600.0)
end

"""
    iau_2006a_crs_cis(date)

Generate the celestial reference system (CRS) to celestial intermediate system
(CIS) transform matrix given the CIP X & Y coordinates and the CIO locator S
for the specified date using the IAU 2006/2000A precession nutation model.

"""
function iau_2006a_crs_cis(date)
    
    pnmat = precession_nutation(tt; model=:iau2006a)
    cio_s = iau_2006_cio_locator(tt, pnmat[3,1:2])
    E = (x*x+y*y) > 0.0 ? atan(y, x) : 0.0
    d = atan(sqrt((x*x+y*y)/(1-(x*x+y*y))))
    Rz(-(E+cio_s))Ry(d)Rz(E)
end

"""
    iau_2006a_crs_trs(date)

Generate the celestial reference system (CRS) to terrestrial reference system
(TRS) transform matrix given given the date, UT1, and polar motion using the
IAU 2006/2000A precession nutation model.

"""
function iau_2006_crs_trs(date, ut1, pole)

    c2i = iau_2006a_crs_cis(date)
    
end

function iau_2006_ecliptic_equator(date, position)
end

function iau_2006_equator_ecliptic(date, position)
end

function iau_2006a_origins(date)
end

"""
    gst_2006a(uct, tt)

Greenwich apparent sidereal time (IAU 2000 and 2006)

"""
function iau_2006a_gst(utc, tt)

    pnmat = precession_nutation(tt; model=:iau2006a)
    cio_s = iau_2006_cio_locator(tt, pnmat[3,1:2])
    iau_2000_era(utc) - eors(pnmat, cio_s)
end

"""
    iau_2006_gmst(ut1, tt)

Greenwich Mean Sidereal Time (GMST) for the specified UT1 and TT dates.
"""
function iau_2006_gmst(date1, date2)

    Δt = (date2 - JD2000)/(100*DAYPERYEAR)
    
    rem2pi(iau_era_2000(date1) + deg2rad(1/3600)*Polynomial(gmst_2006, :Δt)(Δt))
end

function iau_2006a_nutation(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ψ, ϵ = iau_2000a_nutation(date)
    (ψ + ψ*(bj2_2006 + fj2_2006*Δt), ϵ + ϵ*fj2_2006*Δt)
end

function iau_2006_obliquity(date)
    
    Δt = (date - JD2000)/(100*DAYPERYEAR)

    deg2rad(1/3600)*Polynomial(ϵB_2006, :Δt)(Δt)
end

function iau_2006_precession(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    deg2rad(1/3600).*(Polynomial(γB_2006, :Δt)(Δt), Polynomial(ϕB_2006, :Δt)(Δt),
                      Polynomial(ψB_2006, :Δt)(Δt), Polynomial(ϵB_2006, :Δt)(Δt))
end

function iau_2006_tdb_tt(date, ut1, eastlon, u, v)

    Δt = (date - JD2000)/(1000*DAYPERYEAR)
    Δt = ((2448939.5 - JD2000) + 0.123)/(1000*DAYPERYEAR)
    
    # Topocentric terms
    tsol = 2pi*rem(ut1, 1) + eastlon

    # Fundamental arguments (Simon et al. 1994)
    ϵsun  = deg2rad(rem(Polynomial(ϵsun_1994, :Δt)(Δt/3600), 360))
    ϵmsun = deg2rad(rem(Polynomial(ϵmsun_1994, :Δt)(Δt/3600), 360))
    D     = deg2rad(rem(Polynomial(D_1994, :Δt)(Δt/3600), 360))
    ϵju   = deg2rad(rem(Polynomial(ϵju_1994, :Δt)(Δt/3600), 360))
    ϵsa   = deg2rad(rem(Polynomial(ϵsa_1994, :Δt)(Δt/3600), 360))

    wt = sum(topo_1994 .* [
        u * sin(tsol + ϵsun - ϵsa),
        u * sin(tsol - 2*ϵmsun),
        u * sin(tsol - D),
        u * sin(tsol + ϵsun - ϵju),
        u * sin(tsol + 2*ϵsun + ϵmsun),
        v * cos(ϵsun + ϵmsun),
        u * sin(tsol - ϵmsun),
        u * sin(tsol + 2*ϵsun),
        v * cos(ϵsun),
        u * sin(tsol)])
    
    wf = Polynomial([
        sum(tdb_tt_2003_0[1,:] .*
            sin.(tdb_tt_2003_0[3,:] .+ tdb_tt_2003_0[2,:].*Δt))
        sum(tdb_tt_2003_1[1,:] .*
            sin.(tdb_tt_2003_1[3,:] .+ tdb_tt_2003_1[2,:].*Δt))
        sum(tdb_tt_2003_2[1,:] .*
            sin.(tdb_tt_2003_2[3,:] .+ tdb_tt_2003_2[2,:].*Δt))
        sum(tdb_tt_2003_3[1,:] .*
            sin.(tdb_tt_2003_3[3,:] .+ tdb_tt_2003_3[2,:].*Δt))
        sum(tdb_tt_2003_4[1,:] .*
            sin.(tdb_tt_2003_4[3,:] .+ tdb_tt_2003_4[2,:].*Δt))], :Δt)(Δt)
    
    wj = sum(mass_plan_1994_0[1,:] .*
             sin.(mass_plan_1994_0[3,:] .+ mass_plan_1994_0[2,:].*Δt)) +
             mass_plan_1994_2 * Δt^2

    wt + wf + wj
end
