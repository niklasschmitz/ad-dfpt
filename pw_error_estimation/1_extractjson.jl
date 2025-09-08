#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t4 $BN.jl 2>&1 | tee $BN.log
exit $?
=#
using JLD2
using JSON3
using LinearAlgebra
using OrderedCollections
using StaticArrays

workdir = @__DIR__
files = filter(endswith(".jld2"), readdir(joinpath(workdir, "runs")))

function scfres_to_dict(f)
    OrderedDict(
        "lattice"=>f["lattice"],
        "atomic_positions"=>f["atomic_positions"],
        "species"=>f["species"],
        "forces"=>f["forces"],
        "Ecut"=>f["Ecut"],
        "kgrid"=>f["kgrid"],
        "fft_size"=>f["fft_size"],
        "n_electrons"=>f["n_electrons"],
    )
end

for file in files
    if !contains(file, "refined")
        continue
    end

    regex = r"^(.+)_Ecut(\d+)_refined_Ecutlarge([0-9\.]+).jld2$"
    m = match(regex, file)
    if isnothing(m)
        error("File $file does not match expected pattern")
    end

    prefix = m.captures[1]
    Ecutlow = m.captures[2]
    Ecuthigh = m.captures[3]

    data = OrderedDict()

    jldopen(joinpath(workdir, "runs", "$(prefix)_Ecut$(Ecutlow)_scfres.jld2"), "r") do f
        data["base"] = scfres_to_dict(f)
    end
    jldopen(joinpath(workdir, "runs", "$(prefix)_Ecut$(Ecuthigh)_scfres.jld2"), "r") do reference_file
        data["reference"] = scfres_to_dict(reference_file)
        jldopen(joinpath(workdir, "runs", file), "r") do f
            data["refinement"] = OrderedDict(
                "forces_base"=>f["F_base"],
                "forces_refined"=>f["F_refined"],
                "ρ_l1"=>norm(f["ρ"], 1),
                "δρ_l1"=>norm(f["δρ"], 1),
                "ρ_error_l1"=>norm(f["ρ"] - reference_file["ρ"], 1),
                "ρ_plus_δρ_error_l1"=>norm(f["ρ"] + f["δρ"] - reference_file["ρ"], 1),
                "dvol"=>reference_file["dvol"],
            )
        end
    end

    open(joinpath(workdir, "runs", "$prefix.json"), "w") do io
        JSON3.pretty(io, data)
    end
end

