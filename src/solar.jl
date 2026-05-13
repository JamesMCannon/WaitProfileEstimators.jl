"""
    zenithangle(lat_deg, lon_deg, yr, mo, dy, hr, Δτ=nothing) -> Float64
    zenithangle(lat_deg, lon_deg, dt::DateTime, Δτ=nothing)   -> Float64

Solar zenith angle in degrees at geographic latitude `lat_deg`, geographic
longitude `lon_deg` (east positive), and universal time (UT1) given either
as broken-down components `(yr, mo, dy, hr)` with `hr` a fractional hour,
or as a `DateTime` `dt`.

`Δτ` is the difference between terrestrial time and universal time in
seconds. If `nothing`, it is approximated by linear extrapolation; an
error of up to ~30 s has negligible effect on the returned zenith angle.

Maximum error is ≈ 0.0027° (~10 arcsec) for dates between 2010 and 2110.
Refraction is not applied.

This function is *internal* to WaitProfileEstimators and is not exported.

# Reference
Grena, R., 2012. Five new algorithms for the computation of sun position
from 2010 to 2110. *Solar Energy* 86, 1323–1337.
doi: 10.1016/j.solener.2012.01.024
"""
function zenithangle(lat_deg::Real, lon_deg::Real,
                     yr::Integer, mo::Integer, dy::Integer, hr::Real,
                     Δτ = nothing)
    θ, ϕ = deg2rad(lon_deg), deg2rad(lat_deg)

    # ---- Time scale --------------------------------------------------------
    # January and February are treated as months 13 and 14 of the previous
    # year so the integer-month conversion is monotonic.
    if mo <= 2
        mo += 12
        yr -= 1
    end
    t = trunc(Int, 365.25 * (yr - 2000)) +
        trunc(Int,  30.6001 * (mo + 1)) -
        trunc(Int,   0.01 * yr) +
        dy + 0.0416667 * hr - 21958
    if Δτ === nothing
        Δτ = 96.4 + 0.00158 * t
    end
    te = t + 1.1574e-5 * Δτ

    # ---- Heliocentric ecliptic longitude L --------------------------------
    ωa = 0.0172019715
    a  = (3.33024e-2,  3.512e-4,  5.2e-6)
    b  = (-2.0582e-3, -4.07e-5,  -9.0e-7)
    s1, c1 = sincos(ωa * te)
    s2, c2 = 2*s1*c1, (c1 + s1) * (c1 - s1)
    s3, c3 = s2*c1 + c2*s1, c2*c1 - s2*s1
    s = (s1, s2, s3)
    c = (c1, c2, c3)

    β  =  2.92e-5
    dβ = -8.23e-5
    ω = (1.49e-3, 4.31e-3, 1.076e-2, 1.575e-2, 2.152e-2, 3.152e-2, 2.1277e-1)
    d = (1.27e-5, 1.21e-5, 2.33e-5,  3.49e-5,  2.67e-5,  1.28e-5,  3.14e-5)
    φ = (-2.337,  3.065,  -1.533,   -2.358,    0.074,    1.547,   -0.488)

    L = 1.7527901 + 1.7202792159e-2 * te
    for k in 1:3
        L += a[k] * s[k] + b[k] * c[k]
    end
    L += dβ * s1 * sin(β * te)
    for i in 1:7
        L += d[i] * sin(ω[i] * te + φ[i])
    end

    # ---- Nutation correction ----------------------------------------------
    ωn = 9.282e-4
    ν  = ωn * te - 0.8
    sν, cν = sincos(ν)
    Δλ = 8.34e-5 * sν
    λ  = L + π + Δλ
    ϵ  = 4.089567e-1 - 6.19e-9 * te + 4.46e-5 * cν

    # ---- Right ascension / declination ------------------------------------
    sλ, cλ = sincos(λ)
    sϵ, cϵ = sincos(ϵ)
    α = atan(sλ * cϵ, cλ)
    α < 0 && (α += 2π)
    δ = asin(sλ * sϵ)

    # ---- Hour angle --------------------------------------------------------
    H = 1.7528311 + 6.300388099 * t + θ - α + 0.92 * Δλ
    H = mod2pi(H + π) - π

    # ---- Zenith angle with topocentric parallax correction -----------------
    sϕ, cϕ = sincos(ϕ)
    sδ, cδ = sincos(δ)
    cH = cos(H)
    se0 = sϕ * sδ + cϕ * cδ * cH
    ep  = asin(se0) - 4.26e-5 * sqrt(1 - se0^2)
    z   = π/2 - ep
    return rad2deg(z)
end

function zenithangle(lat_deg::Real, lon_deg::Real, dt::DateTime, Δτ = nothing)
    y, m, d = yearmonthday(dt)
    h = hour(dt) + minute(dt)/60 + second(dt)/3600 +
        millisecond(dt)/3_600_000
    return zenithangle(lat_deg, lon_deg, y, m, d, h, Δτ)
end