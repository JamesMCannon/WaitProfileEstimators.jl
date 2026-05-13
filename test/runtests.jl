using WaitProfileEstimators
using Dates
using Test

# Pull in the internal zenithangle for direct testing.
using WaitProfileEstimators: zenithangle

@testset "WaitProfileEstimators.jl" begin

    lola_n = (-115.27,  50.00)   # Northern hemisphere test point
    lola_s = (-115.27, -30.00)   # Southern hemisphere (FIRI seasonal mirror)
    dt_day   = DateTime(2025, 6, 15, 20, 0)   # ≈ local noon at -115° lon
    dt_night = DateTime(2025, 6, 16,  8, 0)   # ≈ local midnight

    @testset "zenithangle (Grena 2012)" begin
        # Sun is up at local noon, down at local midnight.
        @test zenithangle(lola_n[2], lola_n[1], dt_day)   <  90
        @test zenithangle(lola_n[2], lola_n[1], dt_night) >  90

        # June solstice: subsolar point near (23.4°N, sub-solar longitude).
        # At lat = 23.4°N, lon = -180 + 15·(12 - h_UT) ≈ local solar noon,
        # χ should be very small. Pick a UT that lands solar noon at 0° lon.
        @test zenithangle(23.4, 0.0, DateTime(2025, 6, 21, 12, 0)) < 2

        # DateTime and broken-down forms agree exactly.
        @test zenithangle(lola_n[2], lola_n[1], dt_day) ==
              zenithangle(lola_n[2], lola_n[1], 2025, 6, 15, 20.0)

        # Specifying Δτ explicitly produces a finite, similar result.
        z_auto    = zenithangle(lola_n[2], lola_n[1], dt_day)
        z_with_dt = zenithangle(lola_n[2], lola_n[1], dt_day, 69.0)
        @test isfinite(z_with_dt)
        @test abs(z_auto - z_with_dt) < 0.01
    end

    @testset "Return shape and plausibility" begin
        for m in (FIRIFit(), McRaeThomson(), Ferguson())
            (hp, β) = hprime_beta(m, lola_n, dt_day)
            @test isfinite(hp) && isfinite(β)
            @test 50  <  hp <  120     # km
            @test 0.1 <   β <  1.0     # km⁻¹
        end
    end

    @testset "Symbol convenience wrapper" begin
        @test hprime_beta(lola_n, dt_day; method = :firi)     ===
              hprime_beta(FIRIFit(),      lola_n, dt_day)
        @test hprime_beta(lola_n, dt_day; method = :thomson)  ===
              hprime_beta(McRaeThomson(), lola_n, dt_day)
        @test hprime_beta(lola_n, dt_day; method = :ferguson) ===
              hprime_beta(Ferguson(),     lola_n, dt_day)
        @test hprime_beta(lola_n, dt_day; method = :mcraethomson) ===
              hprime_beta(McRaeThomson(), lola_n, dt_day)
        @test_throws ArgumentError hprime_beta(lola_n, dt_day; method = :nope)
    end

    @testset "FIRIFit configuration takes effect" begin
        a = hprime_beta(FIRIFit(),                         lola_n, dt_day)
        b = hprime_beta(FIRIFit(zrange_km = (65.0, 85.0)), lola_n, dt_day)
        @test a != b
    end

    @testset "FIRIFit: southern-hemisphere mirroring runs" begin
        # Just needs to return finite, plausible values — the mirror is an
        # internal convention, not something to assert numerically here.
        hp, β = hprime_beta(FIRIFit(), lola_s, dt_day)
        @test isfinite(hp) && isfinite(β)
        @test 50  <  hp <  120
        @test 0.1 <   β <  1.0
    end

   @testset "McRaeThomson: AM/PM sign of χ" begin
        am = hprime_beta(McRaeThomson(), lola_n, DateTime(2025, 6, 15, 16, 0))
        pm = hprime_beta(McRaeThomson(), lola_n, DateTime(2025, 6, 16,  0, 0))
        @test am != pm
    end

    @testset "McRaeThomson: warns on out-of-range χ" begin
        # Local midnight at lon = -115°, lat = 50° is ~08:00 UTC; χ will be
        # well over 98° there.
        dt_night = DateTime(2025, 6, 16, 8, 0)
        @test_warn r"McRaeThomson"i hprime_beta(McRaeThomson(), lola_n, dt_night)

        # Still returns a (finite) value despite the warning.
        hp, β = hprime_beta(McRaeThomson(), lola_n, dt_night)
        @test isfinite(hp) && isfinite(β)
    end

    @testset "Ferguson: corrected β returned" begin
        # The corrected form should give β ≥ 0.15 at all reasonable inputs
        # (the raw polynomial sits around 0.2–0.5; +0.15 keeps it well above
        # the correction floor).
        for dt in (dt_day, dt_night)
            _, β = hprime_beta(Ferguson(), lola_n, dt)
            @test β > 0.15
        end
    end

    @testset "waitprofile evaluation" begin
        @test waitprofile(70_000, 75.0, 0.45) > 0

        zs  = collect(60_000:1000:90_000)
        nes = waitprofile(zs, 75.0, 0.45)
        @test length(nes) == length(zs)
        @test all(>(0), nes)

        # At z = h', Ne reduces to the prefactor N0 · exp(-0.15 h').
        @test waitprofile(75_000, 75.0, 0.45) ≈ 1.43e13 * exp(-0.15 * 75)
    end
end