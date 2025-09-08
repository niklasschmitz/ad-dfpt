#!/bin/bash

# Collect all .extxyz files in sol58lc
file_list=$(find ../../sol58lc/structures -name "*.extxyz")

# Maximum number of parallel jobs
MAX_JOBS=30

total_jobs=$(echo $file_list | wc -w)
echo "Total jobs: $total_jobs"
echo "Max concurrent jobs: $MAX_JOBS"
echo "-------------------------------------"

# Loop over dataset and enumerate files with index
# and launch jobs in the background
for file in $file_list; do
    # Enforce max concurrent jobs
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 1
    done

    BN=$(basename "$file" .extxyz)
    LOGFILE="$BN.log"
    # Launch job
    (
        echo "Launching $LOGFILE"
        julia --project -t1 0_eval.jl $file > "$LOGFILE" 2>&1
    ) &
done

# Wait for all remaining background jobs to finish
wait

echo "All jobs completed."