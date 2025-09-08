#!/bin/sh
## To be called like `bash 2_linear_pushforward.jl scfres.jld2`.
#=
INPUTFILE=$1
julia --project -t4 $0 $INPUTFILE 2>&1 | tee ${INPUTFILE}_linear_pushforward.log
exit $?
=#
include("lib.jl")
DFTK.setup_threading()
file = ARGS[1]
prefix = basename(file)
save_path = "$(prefix)_linear_pushforward.json"

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

val, grad = @time "implicit differentiation of geometry" begin
    lattice_sensitivity(case, functional, θ_bf; Ecut, kgrid, p0, g_abstol=1e-8,
                        temperature, smearing, kinetic_blowup)
end

# Compute the linearized pushforward as N(a(θ_bf), JΣJᵀ).
# Here J = (∇_θ a)ᵀ is the Jacobian of the map θ ↦ a(θ).
# We have JΣJᵀ = JLLᵀJᵀ = |Lᵀ∇a|^2 = σ²
# Thus we just compute σ := |Lᵀ∇a|
results = Dict()
results["lattice_constant"] = val
results["lattice_constant_grad"] = grad
results["lattice_constant_uncertainty"] = norm(L'grad)
open(save_path, "w") do io
    JSON3.pretty(io, results)
end
