struct Equinox2000
    ic::Vector{Int8}
    fc::Vector{Float64}
end

struct Lunisolar2000
    ic::Vector{Int8}
    #    P::Float64      # (days)
    # the following units are in 0.1 μas
    fc::Vector{Float64}
end

struct Planetary2000
    ic::Vector{Int8}
    #    P::Float64      # (days)
    # the following units are in 0.1 μas
    fc::Vector{Float64}
end

# IAU 2000 Model

include("constants2000.jl")

function iau_2000_equinox_complement(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    longitude = deg2rad(1/3600).*rem.([
        Polynomial(l0_2003A, :Δt)(Δt),
        Polynomial(l1_2003A, :Δt)(Δt),
        Polynomial( F_2003A, :Δt)(Δt),
        Polynomial( D_2000A, :Δt)(Δt),
        Polynomial( Ω_2003A, :Δt)(Δt)], ARCSECPER2PI)
    append!(longitude, [
        Polynomial( lve_2003, :Δt)(Δt),
        Polynomial( lea_2003, :Δt)(Δt),
        Polynomial( lge_2003, :Δt)(Δt)])
    
    sum0, sum1 = 0., 0.
    for term in Iterators.reverse(iau_2000_equinox_0_series)
        sum0 += sum(term.fc.*sincos(sum(term.ic.*longitude)))
    end

    for term in Iterators.reverse(iau_2000_equinox_1_series)
        sum1 += sum(term.fc.*sincos(sum(term.ic.*longitude)))
    end

    deg2rad(1/3600)*(sum0 + sum1*Δt)
end

function iau_2000_gmst(ut, tt)
end

function ephem_position(coef0, coef1, coef2, Δt)

    A0, ϕ0, ν0 = coef0[1,:], coef0[2,:], coef0[3,:]
    A1, ϕ1, ν1 = coef1[1,:], coef1[2,:], coef1[3,:]
    A2, ϕ2, ν2 = coef2[1,:], coef2[2,:], coef2[3,:]

    (sum(A0 .* cos.(ϕ0 .+ ν0 .* Δt)) +
     sum(A1 .* cos.(ϕ1 .+ ν1 .* Δt))*Δt +
     sum(A2 .* cos.(ϕ2 .+ ν2 .* Δt))*Δt^2)
end

function ephem_velocity(coef0, coef1, coef2, Δt)

    A0, ϕ0, ν0 = coef0[1,:], coef0[2,:], coef0[3,:]
    A1, ϕ1, ν1 = coef1[1,:], coef1[2,:], coef1[3,:]
    A2, ϕ2, ν2 = coef2[1,:], coef2[2,:], coef2[3,:]
    
    (-sum(A0 .*  ν0 .* sin.(ϕ0 .+ ν0 .* Δt)) +
      sum(A1 .* (cos.(ϕ1 .+ ν1 .* Δt) .- ν1 .* Δt .* sin.(ϕ1 .+ ν1 .* Δt))) +
      sum(A2 .* (2 .* cos.(ϕ2 .+ ν2 .* Δt) .-
                 ν2 .* Δt .* sin.(ϕ2 .+ ν2 .* Δt)))*Δt)/DAYPERYEAR
end

include("ephemerisDE405.jl")

"""
    iau_2000_earth_position(date; frame=:barycenter)
"""
function iau_2000_position(date; frame=:barycenter)

    @assert frame in [:barycenter, :heliocenter]
    Δt = (date - JD2000)/DAYPERYEAR
    @assert abs(Δt) <= 100.0 "Julian day is not between 1990 and 2100."

    # Sun to Earth ecliptic vector
    position = [
        ephem_position(sun_earth_x_0, sun_earth_x_1, sun_earth_x_2, Δt),
        ephem_position(sun_earth_y_0, sun_earth_y_1, sun_earth_y_2, Δt),
        ephem_position(sun_earth_z_0, sun_earth_z_1, sun_earth_z_2, Δt)]
    velocity = [
        ephem_velocity(sun_earth_x_0, sun_earth_x_1, sun_earth_x_2, Δt),
        ephem_velocity(sun_earth_y_0, sun_earth_y_1, sun_earth_y_2, Δt),
        ephem_velocity(sun_earth_z_0, sun_earth_z_1, sun_earth_z_2, Δt)]

    # Barycenter to Earth ecliptic vector
    if frame == :barycenter
        position .+= [
            ephem_position(bary_sun_x_0, bary_sun_x_1, bary_sun_x_2, Δt),
            ephem_position(bary_sun_y_0, bary_sun_y_1, bary_sun_y_2, Δt),
            ephem_position(bary_sun_z_0, bary_sun_z_1, bary_sun_z_2, Δt)]
        velocity .+= [
            ephem_velocity(bary_sun_x_0, bary_sun_x_1, bary_sun_x_2, Δt),
            ephem_velocity(bary_sun_y_0, bary_sun_y_1, bary_sun_y_2, Δt),
            ephem_velocity(bary_sun_z_0, bary_sun_z_1, bary_sun_z_2, Δt)]
    end
    (iau_2000_bcrs*position, iau_2000_bcrs*velocity)
end

"""

Frame bias and precession, IAU 2000
"""
function iau_2000_precession(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)
    
    ψ = Polynomial(ψ_1977, :Δt)(Δt) + ψ_corr_2000*Δt
    ϵ = Polynomial(ω_1977, :Δt)(Δt) + ϵ_corr_2000*Δt
    χ  = Polynomial(χ_1977, :Δt)(Δt)

    deg2rad(1/3600).*(ψ, ϵ, χ)
end

function iau_2000_tio_locator(date)
end

#    IAU 2000A Model

include("constants2000A.jl")
include("constants2000B.jl")

function iau_2000a_cio_locator(date)
end

function iau_2000a_crs_cis(date)
end

function iau_2000a_crs_trs(date)
end

function iau_2000a_gst(ut, tt)
end

function iau_2000a_nutation(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ###  Luni-solar Nutation
    
    # Fundamental (Delauney) arguments
    longitude = deg2rad(1/3600).*rem.([
        Polynomial(l0_2003A, :Δt)(Δt),
        Polynomial(l1_2000A, :Δt)(Δt),
        Polynomial( F_2003A, :Δt)(Δt),
        Polynomial( D_2000A, :Δt)(Δt),
        Polynomial( Ω_2003A, :Δt)(Δt)], ARCSECPER2PI)

    ψl, ϵl = 0., 0.
    for term in Iterators.reverse(iau_2000a_nutation_lunisolar_series)
        Δr = rem2pi(sum(term.ic.*longitude), RoundToZero)
        ψl += sum((term.fc[1] + term.fc[2]*Δt, term.fc[3]).*sincos(Δr))
        ϵl += sum((term.fc[6], term.fc[4] + term.fc[5]*Δt).*sincos(Δr))
    end

    ###  Planetary Nutation

    longitude = rem2pi.([
        Polynomial(l0_2000A_planet, :Δt)(Δt),
        Polynomial( F_2000A_planet, :Δt)(Δt),
        Polynomial( D_2000A_planet, :Δt)(Δt),
        Polynomial( Ω_2000A_planet, :Δt)(Δt),
        Polynomial(lme_2003, :Δt)(Δt),
        Polynomial(lve_2003, :Δt)(Δt),
        Polynomial(lea_2003, :Δt)(Δt),
        Polynomial(lma_2003, :Δt)(Δt),
        Polynomial(lju_2003, :Δt)(Δt),
        Polynomial(lsa_2003, :Δt)(Δt),
        Polynomial(lur_2003, :Δt)(Δt),
        Polynomial(lne_2003mhb, :Δt)(Δt)], RoundToZero)
    push!(longitude, Polynomial(lge_2003, :Δt)(Δt))

    ψp, ϵp = 0., 0.
    for term in Iterators.reverse(iau_2000a_nutation_planetary_series)
        Δr = rem2pi(sum(term.ic.*longitude), RoundToZero)
        ψp += sum(term.fc[1:2].*sincos(Δr))
        ϵp += sum(term.fc[3:4].*sincos(Δr))
    end

    deg2rad(1e-7/3600).*(ψl + ψp, ϵl + ϵp)
end

function iau_2000a_xys(date)
end

#   IAU 2000B Model

function iau_2000b_cio_locator(date)
end

function iau_2000b_crs_cis(date)
end

function iau_2000b_crs_trs(date)
end

function iau_2000b_gst(ut, tt)
end

function iau_2000b_nutation(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ###  Luni-solar Nutation
    
    # Fundamental (Delauney) arguments from Simon et al. (1994)
    longitude = deg2rad(1/3600).*rem.([
        Polynomial(l0_2000B, :Δt)(Δt),
        Polynomial(l1_2000B, :Δt)(Δt),
        Polynomial( F_2000B, :Δt)(Δt),
        Polynomial( D_2000B, :Δt)(Δt),
        Polynomial( Ω_2000B, :Δt)(Δt)], ARCSECPER2PI)

    ψl, ϵl = 0.0, 0.0
    for term in Iterators.reverse(iau_2000b_nutation_lunisolar_series)
        Δr = rem2pi(sum(term.ic.*longitude), RoundToZero)
        ψl += sum((term.fc[1] + term.fc[2]*Δt, term.fc[3]).*sincos(Δr))
        ϵl += sum((term.fc[6], term.fc[4] + term.fc[5]*Δt).*sincos(Δr))
    end

    deg2rad(1e-7/3600).*(ψl + 1e4*ψ_2000B_planet, ϵl + 1e4*ϵ_2000B_planet)
end

function iau_2000b_xys(date)
end
