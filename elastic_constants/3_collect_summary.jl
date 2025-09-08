using CSV
using DataFrames

function collect_summary(prefix)
    files = filter(endswith(".csv"), readdir(dirname(prefix), join=true))

    # Collect DFPT results
    dfpt_files = filter(startswith(prefix * "_elastic_dfpt"), files)
    dfpt_data = vcat([DataFrame(CSV.File(f)) for f in dfpt_files]...)
    dfpt_data = sort(dfpt_data, :tol, rev=true)
    println(dfpt_data)
    CSV.write(prefix * "_elastic_dfpt_summary.csv", dfpt_data)

    # Collect finite differences results
    fd_files = filter(startswith(prefix * "_elastic_finitediff"), files)
    fd_data = vcat([DataFrame(CSV.File(f)) for f in fd_files]...)
    fd_data = sort(fd_data, [:h, :tol], rev=true)
    println(fd_data)
    CSV.write(prefix * "_elastic_finitediff_summary.csv", fd_data)
end

collect_summary(ARGS[1])
