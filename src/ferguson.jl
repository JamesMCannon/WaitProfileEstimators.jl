"""
    Ferguson()

Ferguson (1980) empirical D-region ionosphere fit. The returned β
includes the +0.15 km⁻¹ correction to the Wait/LongwaveModePropagator
convention. Higher-order terms (sunspot number, magnetic absorption)
are ignored.

# Reference
Ferguson, J.A., 1980. Ionospheric profiles for predicting nighttime
VLF/LF propagation. NOSC/TR-530.
"""
struct Ferguson <: WaitProfileMethod end

function _ferguson(lat_deg::Real, sza_deg::Real, month_num::Integer)
    lat    = deg2rad(lat_deg)
    sza    = deg2rad(sza_deg)
    cos_lat = cos(lat)
    cos_sza = cos(sza)
    month_scalar = cospi(2 * (month_num - 0.5) / 12)
    hprime = 74.37 -
             8.097   * cos_sza +
             5.779   * cos_lat -
             1.213   * month_scalar
    beta   = (0.3849 -
              0.1658  * cos_sza -
              0.08584 * cos_lat +
              0.1296  * month_scalar) + 0.15
    return (hprime, beta)
end

function hprime_beta(::Ferguson, lola::Tuple{<:Real,<:Real}, dt::DateTime)
    lon, lat = lola
    sza = zenithangle(lat, lon, dt)
    return _ferguson(lat, sza, month(dt))
end