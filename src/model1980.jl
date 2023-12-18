struct Nut1980
    ic::Vector{Int8}
    #    P::AbstractFloat      # (days)
    #   The following have units of 0.1 mas.
    fc::Vector{Float64}
end

# Include IAU 1980 luni-solar nutation model constants
include("constants1980.jl")

"""
    iau_1980_nutation(date)

IAU 1980 nutation angles

# Argument
- `date::AbstractFloat`: date of nutation (in TT)

# Returns
- `angles::Tuple{AbstractFloat}`: nutation angles ψ, ϵ (in radians)

The nutation angles are with respect to the ecliptic of date.
"""
function iau_1980_nutation(date::AbstractFloat)

    Δt  = (date - JD2000)/(100*DAYPERYEAR)
    
    l = 2π*rem.([l0_1980t, l1_1980t, F_1980t, D_1980t, Ω_1980t].*Δt, 1.0) .+
        deg2rad.([Polynomial(l0_1980, :Δt)(Δt), Polynomial(l1_1980, :Δt)(Δt),
                  Polynomial( F_1980, :Δt)(Δt), Polynomial( D_1980, :Δt)(Δt),
                  Polynomial( Ω_1980, :Δt)(Δt)]./3600.0)

    ln = vcat([t.n' for t in iau_1980_nutation_series]...)
    la = vcat([t.a' for t in iau_1980_nutation_series]...)
    ϕl = ln*rem2pi.(l, RoundNearest)
    deg2rad.((sum((la[:,1] .+ la[:,2].*Δt).*sin.(ϕl)),
              sum((la[:,3] .+ la[:,4].*Δt).*cos.(ϕl)))./3.6e7)
end

"""
    iau_1980_obliquity(date)

IAU 1980 obliquity model angle

# Arguments
- `date::AbstractFloat`: date of obliquity (in TT)

# Returns
- `angle::AbstractFloat`: angle of obliquity at specified date (in radian)

The angle is the difference between the ecliptic and mean equator of date.
"""
function iau_1980_obliquity(date::AbstractFloat)
    
    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ϵ = deg2rad(1/3600)*Polynomial(ϵ_1980, :Δt)(Δt)
end

"""
    iau_1976_prec(jd1, jd2)

IAU 1976 precession angles

# Arguments
- `date1::Float64`: Starting date of precession (in TDB)
- `date1::Float64`: Ending date of precession (in TDB)

# Returns
- `angles::Tuple{AbstractFloat}`: precession angles ζ, z, θ (in radians)

The rotation matrix is Rz(-z)Ry(θ)Rz(-ζ).

The accumulated precession angles are valid for a limited time span. The
absolute accuracy of the present formulation is <0.1 arcsec for 1960AD -
2040AD, <1 arcsec for 1640AD - 2360AD, and <3 arcsec for 500BC - 3000AD.
The errors are >10 arcsec outside of 1200BC - 3900AD.
"""
function iau_1976_precession(date1::AbstractFloat, date2::AbstractFloat)

    ΔT = (date1 - JD2000)/(100*DAYPERYEAR)
    Δt = (date2 - date1)/(100*DAYPERYEAR)

    ζ = Polynomial([0., Polynomial(ζT_1976, :ΔT)(ΔT),
                     Polynomial(ζA_1976[1:2], :ΔT)(ΔT), ζA_1976[3]], :Δt)(Δt)
    z = Polynomial([0., Polynomial(ζT_1976, :ΔT)(ΔT),
                     Polynomial(zA_1976[1:2], :ΔT)(ΔT), zA_1976[3]], :Δt)(Δt)
    θ = Polynomial([0., Polynomial(θT_1976, :ΔT)(ΔT),
                     Polynomial(θA_1976[1:2], :ΔT)(ΔT), θA_1976[3]], :Δt)(Δt)

    deg2rad(1/3600).*(ζ, z, θ)
end
