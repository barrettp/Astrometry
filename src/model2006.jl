#   IAU 2006 Model

include("constants2006.jl")
include("constants2006F.jl")

"""
    iau_2006_cio_locator(date, coords)

Locate the Celestial Intermediate Origin locator, s, on the equator of the
Celestial Intermediate Pole, given the CIP X and Y coordinates. (IAU 2000A and 2006).
"""
function iau_2006_cio_locator(date, coord)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    longitude = deg2rad(1/3600).*rem.([
        Polynomial(l0_2003A, :Δt)(Δt),
        Polynomial(l1_2003A, :Δt)(Δt),
        Polynomial( F_2003A, :Δt)(Δt),
        Polynomial( D_2003A, :Δt)(Δt),
        Polynomial( Ω_2003A, :Δt)(Δt)], ARCSECPER2PI)
    append!(longitude, [
        Polynomial( lve_2003, :Δt)(Δt),
        Polynomial( lea_2003, :Δt)(Δt),
        Polynomial( lge_2003, :Δt)(Δt)])

    cio_s = cio_s_2006[:]
    for term in Iterators.reverse(iau_2006_equinox_0_series)
        cio_s[1] += sum(term.fc.*sincos(sum(term.ic.*longitude)))
    end

    for term in Iterators.reverse(iau_2006_equinox_1_series)
        cio_s[2] += sum(term.fc.*sincos(sum(term.ic.*longitude)))
    end

    for term in Iterators.reverse(iau_2006_equinox_2_series)
        cio_s[3] += sum(term.fc.*sincos(sum(term.ic.*longitude)))
    end

    for term in Iterators.reverse(iau_2006_equinox_3_series)
        cio_s[4] += sum(term.fc.*sincos(sum(term.ic.*longitude)))
    end

    for term in Iterators.reverse(iau_2006_equinox_4_series)
        cio_s[5] += sum(term.fc.*sincos(sum(term.ic.*longitude)))
    end

    deg2rad(1/3600)*Polynomial(cio_s, :Δt)(Δt) - coord[1]*coord[2]/2.0
end

function iau_2006_cip_xy(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    # lunar, solar, and planetary longitudes
    longitude = deg2rad(1/3600).*rem.([
        Polynomial(l0_2003A, :Δt)(Δt),
        Polynomial(l1_2003A, :Δt)(Δt),
        Polynomial( F_2003A, :Δt)(Δt),
        Polynomial( D_2003A, :Δt)(Δt),
        Polynomial( Ω_2003A, :Δt)(Δt)], ARCSECPER2PI)
    append!(longitude, [
        Polynomial( lme_2003, :Δt)(Δt),
        Polynomial( lve_2003, :Δt)(Δt),
        Polynomial( lea_2003, :Δt)(Δt),
        Polynomial( lma_2003, :Δt)(Δt),
        Polynomial( lju_2003, :Δt)(Δt),
        Polynomial( lsa_2003, :Δt)(Δt),
        Polynomial( lur_2003, :Δt)(Δt),
        Polynomial( lne_2003, :Δt)(Δt),
        Polynomial( lge_2003, :Δt)(Δt)])
    
    # polynomial part of precession-nutation
    xypr = [sum(Polynomial(cip_x_2006, :Δt)(Δt)),
            sum(Polynomial(cip_y_2006, :Δt)(Δt))]

    # !!! The following code can be improved by rearranging the data arrays.

    jaxy = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1] .+ 1
    jasc = [0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0] .+ 1
    japt = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]

    # nutation periodic terms, planetary
    xypl = [0.0, 0.0]
    ialast = length(ampl_2006)
    for ifreq = length(mfapl_2006):-1:1
        sc = sincos(sum(mfapl_2006[ifreq] .* longitude))
        ia = pntr_2006[ifreq + length(mfals_2006)]
        for i = ialast:-1:ia
            xypl[jaxy[i-ia+1]] += ampl_2006[i] * sc[jasc[i-ia+1]] * Δt^japt[i-ia+1]
        end
        ialast = ia-1
    end

    # nutation periodic terms, luni-solar
    xyls = [0.0, 0.0]
    for ifreq = length(mfals_2006):-1:1
        sc = sincos(sum(mfals_2006[ifreq] .* longitude[1:5]))
        ia = pntr_2006[ifreq]
        for i = ialast:-1:ia
            xyls[jaxy[i-ia+1]] += ampl_2006[i] * sc[jasc[i-ia+1]] * Δt^japt[i-ia+1]
        end
        ialast = ia-1
    end

    deg2rad(1/3600).*(xypr .+ 1e-6 .* (xyls .+ xypl))
end

function iau_2006_crs_cis(date)
end

function iau_2006_crs_trs(date)
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
    iau_2000_era(utc) - Eors(pnmat, cio_s)
end

function iau_2006_gmst(utc, tt)
    
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
            sin.(tdb_tt_2003_0[3,:] + tdb_tt_2003_0[2,:].*Δt))
        sum(tdb_tt_2003_1[1,:] .*
            sin.(tdb_tt_2003_1[3,:] + tdb_tt_2003_1[2,:].*Δt))
        sum(tdb_tt_2003_2[1,:] .*
            sin.(tdb_tt_2003_2[3,:] + tdb_tt_2003_2[2,:].*Δt))
        sum(tdb_tt_2003_3[1,:] .*
            sin.(tdb_tt_2003_3[3,:] + tdb_tt_2003_3[2,:].*Δt))
        sum(tdb_tt_2003_4[1,:] .*
            sin.(tdb_tt_2003_4[3,:] + tdb_tt_2003_4[2,:].*Δt))], :Δt)(Δt)
    
    wj = sum(mass_plan_1994_0[1,:] .*
             sin.(mass_plan_1994_0[3,:] .+ mass_plan_1994_0[2,:].*Δt)) +
             mass_plan_1994_2 * Δt^2

    wt + wf + wj
end
