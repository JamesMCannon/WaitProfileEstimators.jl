module WaitProfileEstimators

using Dates
using Statistics
using LinearAlgebra
using FaradayInternationalReferenceIonosphere
import FaradayInternationalReferenceIonosphere as FIRI

# ---- Abstract API ---------------------------------------------------------

"""
    WaitProfileMethod

Abstract supertype for h'/β estimation strategies. Concrete subtypes
provided by this package: [`FIRIFit`](@ref), [`McRaeThomson`](@ref),
[`Ferguson`](@ref).

# Extending

A new estimator is added by defining a subtype and a method of
[`hprime_beta`](@ref):

```julia
struct MyMethod <: WaitProfileMethod end

function WaitProfileEstimators.hprime_beta(::MyMethod,
                                            lola::Tuple{<:Real,<:Real},
                                            dt::DateTime)
    # ... compute h_prime [km] and β [km⁻¹] ...
    return (h_prime, β)
end
```
"""
abstract type WaitProfileMethod end

"""
    hprime_beta(method::WaitProfileMethod, lola, dt) -> (h_prime, β)
    hprime_beta(lola, dt; method::Symbol = :firi) -> (h_prime, β)

Estimate Wait & Spies D-region profile parameters at geographic location
`lola = (lon_deg, lat_deg)` and UTC datetime `dt`.

# Arguments
- `method`: a [`WaitProfileMethod`](@ref) instance — `FIRIFit()`,
  `McRaeThomson()`, or `Ferguson()`. The keyword form accepts the
  symbols `:firi`, `:mcraethomson` (alias `:thomson`), or `:ferguson`,
  each constructed with default settings.
- `lola`: tuple of `(longitude_deg, latitude_deg)`, longitude east-positive.
- `dt`: `DateTime` in UTC (UT1).

# Returns
A `Tuple` `(h_prime, β)` with `h_prime` in km and `β` in km⁻¹, in the
Wait/LongwaveModePropagator convention.

# Examples

```julia
using WaitProfileEstimators, Dates

lola = (-115.27, 50.00)
dt   = DateTime(2025, 6, 15, 20, 0)

hprime_beta(FIRIFit(),                          lola, dt)
hprime_beta(FIRIFit(zrange_km = (65.0, 85.0)),  lola, dt)
hprime_beta(McRaeThomson(),                     lola, dt)
hprime_beta(Ferguson(),                         lola, dt)
hprime_beta(lola, dt; method = :firi)
```
"""
function hprime_beta end

# ---- Internal helpers (not exported) --------------------------------------
include("solar.jl")
include("waitspies.jl")

# ---- Method implementations -----------------------------------------------
include("firi.jl")
include("mcraethomson.jl")
include("ferguson.jl")

# ---- Symbol-form convenience wrapper --------------------------------------

function hprime_beta(lola::Tuple{<:Real,<:Real}, dt::DateTime;
                     method::Symbol = :firi)
    return hprime_beta(_method_from_symbol(method), lola, dt)
end

function _method_from_symbol(s::Symbol)
    s === :firi                              && return FIRIFit()
    (s === :mcraethomson || s === :thomson)  && return McRaeThomson()
    s === :ferguson                          && return Ferguson()
    throw(ArgumentError(
        "Unknown method :$s. Options: :firi, :mcraethomson (alias :thomson), :ferguson."))
end

# ---- Exports --------------------------------------------------------------

export WaitProfileMethod, FIRIFit, McRaeThomson, Ferguson
export hprime_beta
export waitprofile

end # module