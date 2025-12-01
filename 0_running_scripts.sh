#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

# Optional cleanup step
rm -rf /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data/
rm -rf /home/ocdm0351/DPhil/logs/

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"
prev_jobid=""

# Find all scripts starting with a number (1–15, or more), sort numerically, and run them
for script in $(find "$SCRIPT_DIR" -maxdepth 1 -type f -name "*.sh" | sort -V); do
    filename=$(basename "$script")
    num=${filename%%_*}  # extract the leading number, e.g. "1" from "1_script.sh"

    # only run scripts numbered 1–12
    if (( num >= 1 && num <= 12 )); then
        if [ -z "$prev_jobid" ]; then
            jobid=$(sbatch "$script" | awk '{print $4}')
        else
            jobid=$(sbatch --dependency=afterok:$prev_jobid --kill-on-invalid-dep=yes "$script" | awk '{print $4}')
        fi
        echo "Submitted $script as job $jobid"
        prev_jobid=$jobid
    fi
done
