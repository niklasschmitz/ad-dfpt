using AtomsIO
using DFTK
using PseudoPotentialData
using Unitful, UnitfulAtomic
using Plots

pseudopotentials = PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf")

system = load_system(joinpath(@__DIR__, "structures", "GaAs_b3.extxyz"))
model = model_DFT(system; functionals=PBE(), pseudopotentials)

kspacing = 0.15u"Å^-1"
kgrid = kgrid_from_maximal_spacing(model, kspacing)
basis = PlaneWaveBasis(model; Ecut, kgrid)

# Compute an EOS
res = map(-2:2) do i
    new_model = Model(model; lattice=(1 + i * 0.01) * model.lattice)
    new_basis = PlaneWaveBasis(new_model; Ecut, kgrid)
    scfres = self_consistent_field(new_basis; tol=1e-4)
end

plot(
    [(1 + i * 0.01) for i in -2:2],
    [r.energies.total / length(r.basis.model.atoms) for r in res], 
    xlabel="Strain factor",
    ylabel="Energy / atom (Ha)",
    label="DFTK @ dojo.nc.sr.pbe.v0_4_1.standard",
    title="Sol58LC GaAs zincblende PBE",
    m=:o
)
