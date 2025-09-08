using Optim
using ForwardDiff
using LinearAlgebra
using LineSearches
using ComponentArrays

#
# Custom wrapper around Optim for AD support,
# since Optimization.jl is not quite there yet.
#

function optimise(fg!, x0::AbstractVector, θ::AbstractVector; show_trace=false, kwargs...)
    result = Optim.optimize(Optim.only_fg!((E, G, x) -> fg!(E, G, x, θ)),
                            x0, LBFGS(; linesearch=LineSearches.BackTracking()),
                            Optim.Options(; allow_f_increases=true, show_trace, kwargs...))

    if show_trace
        println(result)
    end

    result.minimizer
end

function optimise(fg!, x0::AbstractVector{<:Number}, θ::AbstractVector{T}; kwargs...) where {T <: ForwardDiff.Dual}
    θ_primal = ForwardDiff.value.(θ)
    xstar  = optimise(fg!, x0, θ_primal; kwargs...)
    function wrap_gradient(x, θ)
        G = zeros(promote_type(eltype(x), eltype(θ)), length(x))
        fg!(0.0, G, x, θ)  # TODO can be optimised if the scfres is passed!
        G
    end

    # Let E(x, θ) be our objective (DFT energy) with lattice constant x and θ DFT parameters.
    # Let further S(x, θ) = [∂E/∂x](x, θ) be the derivative (stress). Optim ensures
    # S(x*, θ) = 0, therefore (all taken at (x*, θ)):
    #    0 = dS/dθ = (∂S/∂x) (∂x/∂θ) + ∂S/∂θ
    # or in other words
    #    ∂x/∂θ = - (∂S/∂x)⁻¹ (∂S/∂θ) = - (∂²E/∂x²)⁻¹ (∂²E/∂θ∂x)
    # where we compute the respective Hessian and cross-derivatives using ForwardDiff
    # on the stress function.

    println("Computing hessian ...")
    hess_x = ForwardDiff.jacobian(x -> wrap_gradient(x, θ_primal), xstar)
    hess_x = factorize(hess_x)
    println("Computing cross-derivative ...")
    grads_dual = wrap_gradient(xstar, θ)
    # TODO: try to compute crossderivative with only 1 DFPT solve

    δxstar = ntuple(ForwardDiff.npartials(T)) do α
        - hess_x \ ForwardDiff.partials.(grads_dual, α)
    end

    δxstar = ntuple(ForwardDiff.npartials(T)) do α
        _convert(xstar, δxstar[α])  # Wrap partials back into ComponentVector
    end

    DT = ForwardDiff.Dual{ForwardDiff.tagtype(T)}
    map((xi, δxi...) -> DT(xi, δxi), xstar, δxstar...)
end

_convert(x::AbstractVector, y::AbstractVector) = y
_convert(x::T, y::AbstractVector) where {T <: ComponentVector} = T(y, getaxes(x))
