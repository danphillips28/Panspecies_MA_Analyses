#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

# Optional cleanup step
rm -rf /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data/
rm -rf /home/ocdm0351/DPhil/logs/

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"
prev_jobid=""

# Build a list of scripts sorted ONLY by their numeric prefix
SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"
prev_jobid=""

# Pure numeric sort
sorted_scripts=$(find "$SCRIPT_DIR" -maxdepth 1 -type f -name "*.sh" \
    | awk -F/ '{print $NF}' \
    | awk -F_ '{print $1, $0}' \
    | sort -n \
    | awk '{print $2}')

for script in $sorted_scripts; do
    fullpath="$SCRIPT_DIR/$script"
    num=${script%%_*}

    if (( num >= 1 && num <= 12 )); then
        if [ -z "$prev_jobid" ]; then
            jobid=$(sbatch "$fullpath" | awk '{print $4}')
        else
            jobid=$(sbatch --dependency=afterok:$prev_jobid "$fullpath" | awk '{print $4}')
        fi

        echo "Submitted $script as job $jobid"
        prev_jobid=$jobid

        # FIX: Give SLURM 1 second to register job before next dependency
        sleep 1
    fi
done
