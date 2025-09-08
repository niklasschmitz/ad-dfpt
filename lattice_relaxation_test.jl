using Test
using DftFunctionals
using PseudoPotentialData
using DFTK
using Unitful, UnitfulAtomic
using ComponentArrays
include("lattice_relaxation.jl")
DFTK.setup_threading(; n_blas=1)


function case_Si(p=ComponentArray(; a=10.33))
    name = "Si"
    lattice = p.a / 2 * [[0 1 1.];
                         [1 0 1.];
                         [1 1 0.]]
    Si = ElementPsp(:Si, PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf"))
    atoms     = [Si,        Si        ]
    positions = [ones(3)/8, -ones(3)/8]
    (; name, lattice, atoms, positions, p, Ecut=18, kgrid=(4, 4, 4))
end


@testset "Lattice equilibrium sensitivity w.r.t. PBE params" begin
    θ = parameters(DftFunctional(:gga_x_pbe))
    make_pbe(θ) = [PbeExchange(θ), DftFunctional(:gga_c_pbe)]

    case = case_Si
    Ecut = 18
    kgrid = (2, 2, 2)
    
    popt = @time optimal_lattice(case, make_pbe, θ; Ecut, kgrid)

    v = one.(θ)  # An arbitrary direction in parameter space
    h = 1e-4
    # Note that here we need to be tighter than usual
    aopt1 = optimal_lattice(case, make_pbe, θ - h * v; Ecut, kgrid, p0=popt, g_abstol=1e-9).a
    aopt2 = optimal_lattice(case, make_pbe, θ + h * v; Ecut, kgrid, p0=popt, g_abstol=1e-9).a
    sensitivity_findiff = (aopt2 - aopt1) / 2h

    val, grad = @time lattice_sensitivity(case, make_pbe, θ; Ecut, kgrid, p0=popt)

    @test abs(val - popt.a) < 1e-4
    @test abs(dot(grad, v) - sensitivity_findiff) < h
end

