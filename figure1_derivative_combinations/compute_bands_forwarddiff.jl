using DFTK
using ForwardDiff
using LinearAlgebra


"""This is a wrapper around `compute_bands` which skips perturbation theory of eigenvectors."""
function compute_bands_forwarddiff(basis::PlaneWaveBasis{Float64}, kgrid::DFTK.AbstractKgrid; ρ, kwargs...)
    bands = compute_bands(basis, kgrid; ρ, kwargs...)
    (; bands.basis, bands.eigenvalues)
end

function compute_bands_forwarddiff(basis_dual::PlaneWaveBasis{T}, kgrid::DFTK.AbstractKgrid; ρ, kwargs...) where {T <: ForwardDiff.Dual}
    # Primal pass
    basis_primal = DFTK.construct_value(basis_dual)
    ρ_primal = ForwardDiff.value.(ρ)
    bands_primal = compute_bands(basis_primal, kgrid; ρ=ρ_primal, kwargs...)

    # 1st order perturbation theory of eigenvalues
    basis = PlaneWaveBasis(basis_dual, kgrid)
    ham_dual = energy_hamiltonian(basis, bands_primal.ψ, bands_primal.occupation;
                                  ρ, bands_primal.eigenvalues,
                                  bands_primal.εF).ham
    Hψ_dual = ham_dual * bands_primal.ψ
    δeigenvalues = ntuple(ForwardDiff.npartials(T)) do α
        δHψ = [ForwardDiff.partials.(δHψk, α) for δHψk in Hψ_dual]
        map(bands_primal.ψ, δHψ) do ψk, δHψk
            map(eachcol(ψk), eachcol(δHψk)) do ψnk, δHψnk
                real(dot(ψnk, δHψnk))  # δε_{nk} = <ψnk | δH | ψnk>
            end
        end
    end

    # Convert and combine
    DT = ForwardDiff.Dual{ForwardDiff.tagtype(T)}
    eigenvalues = map(bands_primal.eigenvalues, δeigenvalues...) do εk, δεk...
        map((εnk, δεnk...) -> DT(εnk, δεnk), εk, δεk...)
    end
    
    (; basis, eigenvalues)
end
