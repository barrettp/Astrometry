"""
    equinox(date, model=:iau2006a)
"""

function equinox(date, model=:iau2006a)
    models = [:iau2000a, :iau2000b, :iau2006]
    model_message = "Incorrect model type: options are :iau2000a, :iau2000b, :iau2006a"
    @assert model in models "Incorrect model type"
    
    if model == :iau2000a
        ϵA = iau_1980_obliquity(date) + ϵ_2000_corr*(date - JD2000)/(100*DAYPERYEAR)
        ψn, ϵn = iau_2000a_nutation(date)
        equinox = ψn * cos(ϵA) + equinox_complement(date)
    elseif model == :iau2000b
        ϵA = iau_1980_obliquity(date) + ϵ_2000_corr*(date - JD2000)/(100*DAYPERYEAR)
        ψn, ϵn = iau_2000b_nutation(date)
        equinox = ψn * cos(ϵA) + equinox_complement(date)
    elseif model == :iau2006
        equinox = rem2pi(gst_2006a(0.0, date) - gmst_2006(0.0, date), RoundNearest)
    end

    return equinox
end

"""
    precession_nutation(date; model=:iau2006a)


"""
function precession_nutation(date; model=:iau2006a)
    models = [:iau1980, :iau2000a, :iau2000b, :iau2006a]
    model_message = "Incorrect model type: options are :iau1980, :iau2000a, :iau2000b, :iau2006a"
    @assert model in models "Incorrect model type"
    
    if model == :iau1980
        ζ, z, θ = iau_1976_precession(JD2000, date)
        ϵA = iau_1980_obliquity(date)
        ψ, ϵ = iau_1980_nutation(date)
        rmat = Rx(-(ϵA + ϵ))Rz(-ψ)Rx(ϵA)Rz(-z)Ry(θ)Rz(-ζ)
    elseif model == :iau2000a
        # Frame bias angles: GCRS to 2000.0
        rab, ψb, ϵb, ϵ0 = deg2rad(1/3600).*(icrs_ra_2000, ψ_bias_2000, ϵ_bias_2000, ϵ0_2000)
        ψ, ω, χ = iau_2000_precession(date)
        ϵA = iau_1980_obliquity(date) + deg2rad(1/3600)*ϵ_corr_2000*(date - JD2000)/(100*DAYPERYEAR)
        ψn, ϵn = iau_2000a_nutation(date)
        rmat = Rx(-(ϵA + ϵn))Rz(-ψn)Rx(ϵA) * Rz(χ)Rx(-ω)Rz(-ψ)Rx(ϵ0) * Rx(-ϵb)Ry(ψb*sin(ϵ0))Rz(rab)
    elseif model == :iau2000b
        rab, ψb, ϵb, ϵ0 = deg2rad(1/3600).*(icrs_ra_2000, ψ_bias_2000, ϵ_bias_2000, ϵ0_2000)
        ψ, ω, χ = iau_2000_precession(date)
        ϵA = iau_1980_obliquity(date) + deg2rad(1/3600)*ϵ_corr_2000*(date - JD2000)/(100*DAYPERYEAR)
        ψn, ϵn = iau_2000b_nutation(date)
        rmat = Rx(-(ϵA + ϵn))Rz(-ψn)Rx(ϵA) * Rz(χ)Rx(-ω)Rz(-ψ)Rx(ϵ0) * Rx(-ϵb)Ry(ψb*sin(ϵ0))Rz(rab)
    elseif model == :iau2006a
        γB, ϕB, ψB, ϵB = iau_2006_precession(date)
        ψn, ϵn = iau_2006a_nutation(date)
        rmat = Rx(-(ϵB + ϵn))Rz(-(ψB + ψn))Rx(ϕB)Rz(γB)
    end

    return rmat
end
