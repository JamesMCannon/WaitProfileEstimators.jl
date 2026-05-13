"""
    FIRIFit(; zrange_km=(60.0, 90.0), f10_7=(75, 200))

Estimation method that samples a FIRI (Faraday International Reference
Ionosphere) electron density profile at the given location/time and fits
Wait & Spies (h′, β) by OLS on `log(Ne)` vs. altitude over `zrange_km`.

Southern-hemisphere inputs are handled by mirroring latitude and shifting
month by six (per Friedrich et al. 2018), since FIRI is northern-only.

# Keyword arguments
- `zrange_km`: altitude window (km) for the OLS fit. Default `(60.0, 90.0)`.
- `f10_7`: F10.7 solar flux range passed to FIRI. Default `(75, 200)`.
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