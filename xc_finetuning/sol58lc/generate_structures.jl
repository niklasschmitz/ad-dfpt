using CSV
using DataFrames
using ASEconvert


const strukturbericht2ase = Dict(
    "fcc"=>"fcc",
    "bcc"=>"bcc",
    "diamond"=>"diamond",
    "b1"=>"rocksalt",
    "b2"=>"cesiumchloride",
    "b3"=>"zincblende",
)

df = DataFrame(CSV.File(joinpath(@__DIR__, "Sol58LC.csv")))

for row in eachrow(df)
    sysname = String(row["sysnames"])
    a0_exp = row[" ref"]
    a0_pbe = row[" PBE"] + a0_exp  # The original file gives differences to expt
    name, crystalstructure = split(sysname, "_")

    # Convert to Strukturbericht notation to ASE convention
    crystalstructure = strukturbericht2ase[crystalstructure]

    # Build bulk with ASE
    system = ase.build.bulk(name, crystalstructure; a=a0_pbe)
    system.info["a0_exp"] = a0_exp
    system.info["a0_pbe"] = a0_pbe
    system.info["name"] = name
    system.info["crystalstructure"] = crystalstructure
    ase.io.write(joinpath(@__DIR__, "structures", "$sysname.extxyz"), system)
end
