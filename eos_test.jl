using Test
using ForwardDiff
include("eos.jl")


@testset "Birch-Murnaghan EOS of Li-BCC" begin
    # All-electron wien2k results for Li-BCC
    # Raw data from https://archive.materialscloud.org/record/2023.81
    # Extracted using https://github.com/niklasschmitz/pseudopotentials-uncertainty/blob/df4f74d1ebbec2d8694f44d26c4b80bdc3310b6d/external-data/get_energy_volume_data_from_aiida.py
    volumes = [
        19.049974180134118,
        19.455287291015523,
        19.860606964155465,
        20.26592530958009 ,
        20.67123500533444 ,
        21.076557453113914,
        21.481885520814668,
    ]  # in Å^3
    energies = [
        -204.68497361791,
        -204.68694304199,
        -204.68807476354,
        -204.68843926006,
        -204.68810061436,
        -204.68711719486,
        -204.68554206377,
    ]  # in eV

    # Computed using https://github.com/aiidateam/acwf-verification-scripts/blob/main/3-analyze/eos_utils/eosfit_31_adapted.py
    target = (;
        volume0=20.267495555897856,
        e0=-204.68843936504592,
        bulk_modulus0=0.08668095123425397,
        bulk_deriv0=3.3468965122352063,
    )

    for root_method in [Polynomials.roots, roots_abc]
        params = eos_birch_murnaghan_fit(volumes, energies; root_method)

        @test isapprox(params.volume0, target.volume0)
        @test isapprox(params.e0, target.e0)
        @test isapprox(params.bulk_modulus0, target.bulk_modulus0)
        @test isapprox(params.bulk_deriv0, target.bulk_deriv0)
    end

    @testset "ForwardDiff through bm_fit with roots_abc" begin
        using ForwardDiff
        f(x) = eos_birch_murnaghan_fit(volumes, x; root_method=roots_abc).bulk_modulus0
        g = ForwardDiff.gradient(f, energies)

        h = 1e-4
        g2 = (f(energies + h/2 * energies) - f(energies -  h/2 * energies)) / h
        @test isapprox(g'energies, g2, atol=1e-8)
    end
end
