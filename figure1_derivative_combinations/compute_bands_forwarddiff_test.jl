using DFTK
using PseudoPotentialData
using ForwardDiff
using AtomsBuilder
using DifferentiationInterface
include("compute_bands_forwarddiff.jl")
using Test
using Brillouin
using Unitful, UnitfulAtomic

@testset "compute_bands_forwarddiff against finite differences" begin
    pseudopotentials = PseudoFamily("cp2k.nc.sr.lda.v0_1.semicore.gth")
    model = model_DFT(bulk(:Si); pseudopotentials, functionals=LDA(), kinetic_blowup=BlowupCHV())
    kpath = irrfbz_path(model)
    kinter = Brillouin.interpolate(kpath, density=austrip(40u"bohr"))
    kgrid_bs = ExplicitKpoints(kinter)

    function f(h)
        lattice = (1 + h) * model.lattice
        new_model = Model(model; lattice)
        new_basis = PlaneWaveBasis(new_model; Ecut=15, kgrid=[4, 4, 4])
        scfres = self_consistent_field(new_basis, tol=1e-8)
        band_data = compute_bands_forwarddiff(new_basis, kgrid_bs; scfres.ρ, tol=1e-8)
        band_data.eigenvalues[1]
    end

    y0 = f(0.)
    y1, dy1 = value_and_derivative(f, AutoForwardDiff(), 0.)
    dy2 = let h = 1e-5
        (f(h) - f(-h)) / 2h 
    end

    @test maximum(abs, y0 - y1) < 1e-9
    @test maximum(abs, dy1 - dy2) < 1e-5
end
