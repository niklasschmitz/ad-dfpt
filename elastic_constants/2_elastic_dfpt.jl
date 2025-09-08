#!/bin/sh
## To be called like `bash 2_elastic_dfpt.jl scfres.jld2 1e-10`.
#=
INPUTFILE=$1
TOL=$2
julia --project -t4 $0 $INPUTFILE $TOL 2>&1 | tee ${INPUTFILE}_elastic_dfpt_tol${TOL}.log
exit $?
=#
include("lib.jl")
using DataFrames
using CSV
DFTK.setup_threading()
file = ARGS[1]
tol_str = ARGS[2]
tol = parse(Float64, tol_str)
prefix = file

scfres0 = load_scfres(file)
model0 = scfres0.basis.model
kgrid  = scfres0.basis.kgrid
Ecut   = scfres0.basis.Ecut
v0 = model0.unit_cell_volume

@info now() file v0 Ecut kgrid

function symmetries_from_strain(model0, voigt_strain)
    lattice = DFTK.voigt_strain_to_full(voigt_strain) * model0.lattice
    model = Model(model0; lattice, symmetries=true)
    model.symmetries
end

function stress_from_strain(voigt_strain; symmetries, Ecut, kgrid, tol)
    lattice = DFTK.voigt_strain_to_full(voigt_strain) * model0.lattice
    model = Model(model0; lattice, symmetries)
    basis = PlaneWaveBasis(model; Ecut, kgrid)
    scfres = self_consistent_field(basis; tol)
    DFTK.full_stress_to_voigt(compute_stresses_cart(scfres))
end

strain0 = zeros(6)  # The base point is zero strain relative to the relaxed structure
strain_pattern = [1., 0., 0., 1., 0., 0.]  # should yield [c11, c12, c12, c44, 0, 0]

# For elastic constants beyond the bulk modulus, symmetry-breaking strains
# are required. That is, the symmetry group of the crystal is reduced.
# Here we simply precompute the relevant subgroup by applying the automatic
# symmetry detection (spglib) to a finitely perturbed crystal.
symmetries_strain = symmetries_from_strain(model0, 0.01 * strain_pattern)

@info now() length(model0.symmetries) length(symmetries_strain) strain_pattern

f(voigt_strain) = stress_from_strain(voigt_strain; symmetries=symmetries_strain, Ecut, kgrid, tol)
stress, (dstress,) = value_and_pushforward(f, AutoForwardDiff(), strain0, (strain_pattern,))

c11 = ustrip(uconvert(u"GPa", dstress[1] * u"hartree" / u"bohr"^3))
c12 = ustrip(uconvert(u"GPa", dstress[2] * u"hartree" / u"bohr"^3))
c44 = ustrip(uconvert(u"GPa", dstress[4] * u"hartree" / u"bohr"^3))

results = (; tol, Ecut, c11, c12, c44, stress, dstress)
@info now() results...
df = DataFrame()
push!(df, results)
CSV.write(prefix * "_elastic_dfpt_tol$(tol_str).csv", df)
@info "Finished." now()
