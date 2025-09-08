#!/bin/sh
#=
julia -t1 --project -e "include(\"$0\"); main()"
exit $?
=#
using JLD2
using JSON3
using DFTK

function main()
    # Go through all jld2 and extract data for plotting

    files_to_process = filter(endswith(".jld2"), readdir())
    !isempty(files_to_process) && println("Processing files: ")
    for file in files_to_process
        println("   - ", file)
        prefix, _ = splitext(file)
        if endswith(file, "_scfres.jld2")
            scfres = load_scfres(file)
            y = load(file, "y")
            θ = load(file, "θ")
            bands = load(file, "bands")
            extra_data = Dict(
                "y" => Dict(
                    "energies" => y.energies,
                    "forces" => y.forces,
                    "occupation" => y.occupation,
                    "εF" => y.εF,
                    "eigenvalues" => y.eigenvalues,
                    "eigenvalues_bs" => y.eigenvalues_bs,
                    "ρ" => y.ρ[:, :, 1, 1],  # Take only 2D slice along z=0 plane
                    "r" => r_vectors_cart(scfres.basis)[:, :, 1, 1],
                ),
                "θ" => θ,
                "bands" => bands,
            )
            save_scfres(prefix * ".json", scfres; extra_data, save_ρ=false)
        else
            (jac, θ) = load(file, "jac", "θ")
            data = Dict(
                "δenergies" => jac.energies,
                "δforces" => jac.forces,
                "δoccupation" => jac.occupation,
                "δεF" => jac.εF,
                "δeigenvalues" => jac.eigenvalues,
                "δeigenvalues_bs" => jac.eigenvalues_bs,
                "δρ" => jac.ρ[:, :, 1, 1],  # Take only 2D slice along z=0 plane
                "θ" => θ,
            )
            open(prefix * "_plot.json", "w") do fp
                JSON3.write(fp, data)
            end
        end
    end
    nothing
end
