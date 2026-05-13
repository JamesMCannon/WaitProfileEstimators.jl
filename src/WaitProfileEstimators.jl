module WaitProfileEstimators

using Dates
using Statistics
using LinearAlgebra
using FaradayInternationalReferenceIonosphere
import FaradayInternationalReferenceIonosphere as FIRI

# ---- Abstract API ---------------------------------------------------------

"""
    WaitProfileMethod

Abstract supertype for h′/β estimation strategies. Subtypes provided here:
[`FIRIFit`](@ref), [`McRaeThomson`](@ref), [`Ferguson`](@ref).

Extending the package with a new method is just two steps:

    struct MyMethod <: WaitProfileMethod end
    WaitProfileEstimators.hprime_beta(::MyMethod, lola, dt) = ...
"""
abstract type WaitProfileMethod end

"""
    hprime_beta(method::WaitProfileMethod, lola, dt) -> (h_prime, β)
    hprime_beta(lola, dt; method::Symbol = :firi) -> (h_prime, β)

Estimate Wait & Spies D-region profile parameters at the geographic
location `lola = (lon_deg, lat_deg)` and UTC datetime `dt`. Returns a
`Tuple` `(h_prime, β)` with `h_prime` in km and `β` in km⁻¹, in the
Wait/LongwaveModePropagator convention.

The first form dispatches on the concrete method type (recommended; allows
configuration via the method's constructor). The second is a convenience
wrapper accepting `:firi`, `:mcraethomson` (alias `:thomson`), or
`:ferguson`, using default settings for each.

# Examples

```julia
using WaitProfileEstimators, Dates

lola = (-115.27, 50.00)
dt   = DateTime(2025, 6, 15, 20, 0)

hprime_beta(FIRIFit(), lola, dt)
hprime_beta(McRaeThomson(), lola, dt)
hprime_beta(Ferguson(), lola, dt)

hprime_beta(lola, dt; method = :firi)
hprime_beta(FIRIFit(zrange_km = (65.0, 85.0)), lola, dt)
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