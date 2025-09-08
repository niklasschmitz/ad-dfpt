using AtomsIO
using DFTK
using PseudoPotentialData
using Dates
using JSON3
using JLD2
using Random
using LinearAlgebra
using ForwardDiff
include(joinpath(@__DIR__, "../xc_finetuning/sol58lc/load.jl"))
include(joinpath(@__DIR__, "../lattice_relaxation.jl"))
include(joinpath(@__DIR__, "../beef.jl"))
include(joinpath(@__DIR__, "../beefvdw.jl"))
include(joinpath(@__DIR__, "../eos.jl"))


"""Helper function for ForwardDiff for 1st-order Hellmann-Feynman derivatives."""
function total_energy(basis::PlaneWaveBasis; kwargs...)
    scfres = self_consistent_field(basis; kwargs...)
    scfres.energies.total
end

"""ForwardDiff rule for 1st-order Hellmann-Feynman total energy derivatives."""
function total_energy(basis::PlaneWaveBasis{T}; kwargs...) where {T <: ForwardDiff.Dual}
    basis_primal = DFTK.construct_value(basis)
    scfres = self_consistent_field(basis_primal; kwargs...)
    ρ = compute_density(basis, scfres.ψ, scfres.occupation)
    (; energies) = DFTK.energy(basis, scfres.ψ, scfres.occupation; ρ, scfres.eigenvalues, scfres.εF)
    energies.total
end

function replace_xc(model::Model, xc::Xc)
    terms = map(model.term_types) do term
        term isa Xc ? xc : term
    end
    Model(model; terms)
end
