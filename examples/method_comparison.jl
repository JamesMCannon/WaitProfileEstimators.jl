# Generates assets/method_comparison.png — the comparison plot embedded
# in README.md. Run from the package root:
#
#     julia --project=examples examples/method_comparison.jl
#
# Requires Plots and FaradayInternationalReferenceIonosphere in the
# `examples/` environment (see examples/Project.toml).

using WaitProfileEstimators
using FaradayInternationalReferenceIonosphere
import FaradayInternationalReferenceIonosphere as FIRI
using Dates
using Plots

lola = (-115.27, 50.00)
dt   = DateTime(2025, 6, 15, 20, 0)

chi     = WaitProfileEstimators.zenithangle(lola[2], lola[1], dt)
ne_firi = firi(chi, lola[2]; month = month(dt))
z_km    = FIRI.ALTITUDE ./ 1000

methods_table = [
    ("FIRI fit",      FIRIFit()),
    ("McRae–Thomson", McRaeThomson()),
    ("Ferguson",      Ferguson()),
]

plt = plot(xaxis = :log, xlims = (1e6, 1e12), ylims = (60, 95),
           xlabel = "Electron density (m⁻³)",
           ylabel = "Altitude (km)",
           title  = "Wait & Spies estimators at $(lola), $(dt) UTC",
           legend = :bottomright,
           dpi    = 150,
           size   = (700, 500))

plot!(plt, ne_firi, z_km; label = "FIRI sample", lw = 2, color = :goldenrod)

zgrid = 60_000:100:95_000
for (name, m) in methods_table
    hp, β = hprime_beta(m, lola, dt)
    plot!(plt, waitprofile(zgrid, hp, β), zgrid ./ 1000;
          label = "$name  (h' = $(round(hp; digits=2)), β = $(round(β; digits=3)))",
          lw = 2)
end

mkpath(joinpath(@__DIR__, "..", "assets"))
savefig(plt, joinpath(@__DIR__, "..", "assets", "method_comparison.png"))
println("Wrote assets/method_comparison.png")