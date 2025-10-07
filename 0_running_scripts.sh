#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

# Seems to be helpful to remove this before starting: Don't know why exactly..
rm -r /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data/
rm -r /home/ocdm0351/DPhil/logs/

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"
jobid=""

prev_jobid=""
for script in $(ls [1-5]*.sh | sort -n); do
    if [ -z "$prev_jobid" ]; then
        # First job has no dependency
        jobid=$(sbatch "$script" | awk '{print $4}')
    else
        # Later jobs depend on the previous one finishing
        jobid=$(sbatch --dependency=afterok:$prev_jobid --kill-on-invalid-dep=yes "$script" | awk '{print $4}')
    fi
    echo "Submitted $script as job $jobid"
    prev_jobid=$jobid
done
