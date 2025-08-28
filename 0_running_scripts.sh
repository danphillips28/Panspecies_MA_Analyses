#!/bin/bash

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"
jobid=""

prev_jobid=""
for script in $(ls [1-5]*.sh | sort -n); do
    if [ -z "$prev_jobid" ]; then
        # First job has no dependency
        jobid=$(sbatch "$script" | awk '{print $4}')
    else
        # Later jobs depend on the previous one finishing
        jobid=$(sbatch --dependency=afterok:$prev_jobid "$script" | awk '{print $4}')
    fi
    echo "Submitted $script as job $jobid"
    prev_jobid=$jobid
done
