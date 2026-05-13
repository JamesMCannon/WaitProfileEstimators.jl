"""
    FIRIFit(; zrange_km = (60.0, 90.0), f10_7 = (75, 200))

Estimation method that samples a FIRI (Faraday International Reference
Ionosphere) electron density profile at the given location and time and
fits Wait & Spies `(h', β)` by ordinary least squares (OLS) on `log(Ne)` vs.
altitude over `zrange_km`.

Southern-hemisphere inputs are handled by mirroring latitude across the
equator and shifting the month by six, per Friedrich et al. (2018).

# Keyword arguments
- `zrange_km`: altitude window (km) for the OLS fit. Default `(60.0, 90.0)`.
- `f10_7`: F10.7 solar flux range passed to FIRI. Default `(75, 200)`.

# Reference
Friedrich, M., Pock, C., & Torkar, K. (2018). FIRI-2018, an updated
empirical model of the lower ionosphere. *J. Geophys. Res. Space Physics*
123, 6737–6751.

See also: [`hprime_beta`](@ref), [`McRaeThomson`](@ref), [`Ferguson`](@ref).
"""
Base.@kwdef struct FIRIFit <: WaitProfileMethod
    zrange_km::Tuple{Float64,Float64} = (60.0, 90.0)
    f10_7::Tuple{Int,Int}             = (75, 200)
end

function _firi_profile_at(lola::Tuple{<:Real,<:Real}, dt::DateTime, f10_7)
    lon, lat = lola
    chi = zenithangle(lat, lon, dt)

    eff_lat, eff_mo = lat, month(dt)
    if lat < 0
        eff_lat = -lat
        eff_mo  = mod1(eff_mo + 6, 12)
    end
    return firi(chi, eff_lat; f10_7 = f10_7, month = eff_mo)
end

function hprime_beta(m::FIRIFit, lola::Tuple{<:Real,<:Real}, dt::DateTime)
    z  = FIRI.ALTITUDE
    ne = _firi_profile_at(lola, dt, m.f10_7)
    f  = _fit_waitspies_ols(z, ne; zrange_km = m.zrange_km)
    return (f.h_prime, f.β)
end