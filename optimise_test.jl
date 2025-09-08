using Test
using DifferentiationInterface
using LinearAlgebra
include("optimise.jl")

@testset "Implicit Differentiation with ForwardDiff" begin
    f(x, θ) = sum(abs2, x - θ)

    function fg!(_, G, x, θ)
        isnothing(G) && return f(x, θ)
        T = promote_type(eltype(x), eltype(θ))
        (E, _) = value_and_gradient!(x -> f(x, θ), G, AutoForwardDiff(), T.(x))
        E
    end

    N = 10
    θ = ones(N)
    x0 = zeros(N)

    implicit_function(θ) = optimise(fg!, x0, θ)

    x1 = implicit_function(θ)
    x2, jac = value_and_jacobian(implicit_function, AutoForwardDiff(), θ)
    @show x1
    @show x2
    @assert norm(x1 - x2) < 1e-14
    @assert jac ≈ I(N)
end

