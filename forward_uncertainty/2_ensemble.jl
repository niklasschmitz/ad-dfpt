#!/bin/sh
## To be called like `bash 2_ensemble.jl scfres.jld2`.
#=
INPUTFILE=$1
julia --project -t4 $0 $INPUTFILE 2>&1 | tee ${INPUTFILE}_ensemble.log
exit $?
=#
include("lib.jl")
DFTK.setup_threading()
file = ARGS[1]
prefix = basename(file)
save_path = "$(prefix)_ensemble.json"

# Load pre-relaxed system
scfres0 = load_scfres(file)
model0 = scfres0.basis.model

lattice0  = model0.lattice
atoms     = model0.atoms
positions = model0.positions
a0 = load(file, "a_opt")
p0 = ComponentVector(; a=a0)

temperature    = model0.temperature
smearing       = model0.smearing
kinetic_blowup = model0.term_types[1].blowup 

kgrid = scfres0.basis.kgrid
Ecut  = scfres0.basis.Ecut

function case(p)
    lattice = (p.a / a0) * lattice0
    (; lattice, atoms, positions)
end

# Load functional and parameter distribution θ ~ N(θ_bf, Σ) where Σ = LLᵀ
functional = make_beef
(; θ_bf, L) = params_beef_2005()

# Sample ensemble
num_samples = 10
rng = MersenneTwister(1234)
α_ens = randn(rng, length(θ_bf), num_samples)
θ_ens = [θ_bf + L * αi  for αi in eachcol(α_ens)]

results = Dict()
results["theta_ens"] = θ_ens

# Baseline 1: Ensemble with NSCF EOS
volume_scaling_range = [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06]

# Compute self-consistent results across volumes for mean parameter only
model_bf = replace_xc(model0, Xc(functional(θ_bf)))
scfres_curve = map(volume_scaling_range) do scale
    model = Model(model_bf; lattice=cbrt(scale) * lattice0)
    basis = PlaneWaveBasis(model; Ecut, kgrid)
    scfres = self_consistent_field(basis; tol=1e-9)
end

results["ensemble_nscf"] = Dict()
results["ensemble_nscf"]["volume_scaling_range"] = volume_scaling_range
results["ensemble_nscf"]["meanparameter_total_energies"] = [scfres.energies.total / length(scfres.basis.model.atoms) for scfres in scfres_curve]
results["ensemble_nscf"]["total_energies"] = []
results["ensemble_nscf"]["eos_fits"] = []
results["ensemble_nscf"]["lattice_constants"] = []
@time "Ensemble with NSCF EOS" begin
    for (i, θi) in enumerate(θ_ens)
        @show i
        # Compute approximate total energies non-self-consistently
        # TODO: this can be made faster if we only instantiate the TermXc in isolation
        model_i = replace_xc(model_bf, Xc(functional(θi)))
        nscf_energies = map(scfres_curve) do scfres
            model = Model(model_i; lattice=scfres.basis.model.lattice)
            basis = PlaneWaveBasis(model; Ecut, kgrid)
            (; energies) = DFTK.energy(basis, scfres.ψ, scfres.occupation;
                                       scfres.ρ, scfres.eigenvalues, scfres.εF)
            energies.total / length(model.atoms)
        end

        # Fit Equation Of State
        eos = eos_birch_murnaghan_fit(volume_scaling_range, nscf_energies)

        # Convert equilibrium volume (relative to mean-parameter equilibrium volume)
        # to equilibrium lattice constant prediction.
        # NOTE: The scaling range is with respect to the volume of
        #       the relaxed a0 of θ_bf (not of theta_i).
        a = cbrt(eos.volume0) * a0

        push!(results["ensemble_nscf"]["total_energies"], nscf_energies)
        push!(results["ensemble_nscf"]["eos_fits"], eos)
        push!(results["ensemble_nscf"]["lattice_constants"], a)
        open(save_path, "w") do io
            JSON3.pretty(io, results)
        end
    end
end

# Baseline 2: Ensemble with full relaxations
results["ensemble"] = Dict()
results["ensemble"]["lattice_constants"] = []
@time "Ensemble with full relaxations" begin
    for (i, θi) in enumerate(θ_ens)
        @show i
        a = optimal_lattice(case, functional, θi; Ecut, kgrid, p0, g_abstol=1e-8,
                            temperature, smearing, kinetic_blowup).a
        push!(results["ensemble"]["lattice_constants"], a)
        open(save_path, "w") do io
            JSON3.pretty(io, results)
        end
    end
end


