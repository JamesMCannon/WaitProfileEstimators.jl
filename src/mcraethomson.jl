"""
    McRaeThomson()

Daytime midlatitude Wait & Spies parameters `(h', β)` from the
McRae & Thomson (2000) polynomial fits in *signed* solar zenith angle
(negative AM, positive PM). The sign of χ is recovered by
finite-differencing the zenith angle across a one-minute interval.

# Validity
Calibrated near solar minimum for `|χ| ≲ 75°` (≈ 1.71 rad). Inputs
outside this range cause the polynomials to extrapolate and produce
values that should not be trusted; night is not represented. Calls with
`|χ| > 75°` emit a one-time warning per call site, but a value is still
returned — the caller is responsible for restricting use to the valid
range.

# Reference
McRae, W. M., & Thomson, N. R. (2000). VLF phase and amplitude: daytime
ionospheric parameters. *J. Atmos. Sol.-Terr. Phys.* 62, 609–618.

See also: [`hprime_beta`](@ref), [`FIRIFit`](@ref), [`Ferguson`](@ref).
"""
struct McRaeThomson <: WaitProfileMethod end

const MCRAETHOMSON_CHI_LIMIT_DEG = 75.0

function _mcraethomson_poly(chi_rad::Real)
    X = chi_rad
    h_prime = evalpoly(X, (70.55, 0.045, 1.270, -0.774,
                            3.115, 0.595, -1.491, -0.122, 0.268))
    β       = evalpoly(X, (0.395, 0.005, -0.043, 0.0075,
                          -0.0191, -0.0054))
    return (h_prime, β)
end

function hprime_beta(::McRaeThomson, lola::Tuple{<:Real,<:Real}, dt::DateTime)
    lon, lat = lola
    chi_now   = zenithangle(lat, lon, dt)
    chi_after = zenithangle(lat, lon, dt + Minute(1))
    s = chi_after >= chi_now ? +1 : -1
    if chi_now > MCRAETHOMSON_CHI_LIMIT_DEG
        @warn """
              McRaeThomson called with |χ| = $(round(chi_now; digits=2))° > \
              $(MCRAETHOMSON_CHI_LIMIT_DEG)°. The McRae & Thomson (2000) \
              fits are valid for daytime conditions only; the returned \
              (h', β) is an extrapolation and should not be trusted.
              """ maxlog=1 _id=:wpe_mcraethomson_chi_oor
    end
    return _mcraethomson_poly(s * deg2rad(chi_now))
end