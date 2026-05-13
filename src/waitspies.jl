# Wait & Spies constants
const N0_WAIT = 1.43e13      # m⁻³, prefactor in waitprofile
const NU0_KM  = 0.15         # km⁻¹, collision-frequency scale

"""
    waitprofile(z_m, h_prime_km, β_per_km) -> Real
    waitprofile(z_m::AbstractVector, h_prime_km, β_per_km) -> Vector

Evaluate the Wait & Spies two-parameter D-region electron density profile:

    Ne(z) = 1.43e13 · exp(-0.15 h′) · exp((β - 0.15)(z_km - h′))   [m⁻³]

where `z_m` is altitude in **metres**, `h_prime_km` is the reference
height in km, and `β_per_km` is the sharpness in km⁻¹.
"""
function waitprofile(z_m::Real, h_prime_km::Real, β_per_km::Real)
    zkm = z_m / 1000
    return N0_WAIT * exp(-NU0_KM * h_prime_km) *
                     exp((β_per_km - NU0_KM) * (zkm - h_prime_km))
end

function waitprofile(z_m::AbstractVector, h_prime_km::Real, β_per_km::Real)
    return waitprofile.(z_m, h_prime_km, β_per_km)
end

"""
    _fit_waitspies_ols(z_m, ne; zrange_km=(60.0, 90.0)) -> NamedTuple

Internal: OLS fit of Wait & Spies (h′, β) to a discrete (z, Ne) profile,
in log(Ne) vs. altitude (km) over the window `zrange_km`. Returns
`(h_prime, β, log_rms, n)` for diagnostic purposes; the public
`hprime_beta` interface drops the diagnostics and returns just `(h′, β)`.
"""
function _fit_waitspies_ols(z_m::AbstractVector, ne::AbstractVector;
                            zrange_km::Tuple{<:Real,<:Real} = (60.0, 90.0))
    length(z_m) == length(ne) ||
        throw(DimensionMismatch("z_m and ne lengths differ"))
    zkm  = z_m ./ 1000
    mask = (zkm .>= zrange_km[1]) .& (zkm .<= zrange_km[2]) .&
           isfinite.(ne) .& (ne .> 0)
    sum(mask) >= 3 ||
        error("Fewer than 3 valid points in zrange_km = $(zrange_km)")

    x = zkm[mask]
    y = log.(ne[mask])
    X = hcat(ones(length(x)), x)
    b, a = (X \ y)

    β       = a + NU0_KM
    h_prime = (log(N0_WAIT) - b) / β
    rms     = sqrt(mean((y .- (b .+ a .* x)).^2))
    return (h_prime = h_prime, β = β, log_rms = rms, n = length(x))
end