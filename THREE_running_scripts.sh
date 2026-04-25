#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"

# List the scripts explicitly, in the order you want them submitted
scripts=(
"17_Associating_TimeSeries_Expression_with_Mutation_Rate_Submitter.sh"
"18_Associating_Expression_with_Mutation_Rate_Submitter.sh"
"19_Associating_Epigenetics_and_Expression_etc_with_Mutation_Rate_Submitter.sh"
"19_Calculating_CelegansPpacificus_dNdS_Submitter.sh"
"20_Celegans_to_Ppacificus_Expression_Changing_Submitter.sh")

prev_jobid=""

for script in "${scripts[@]}"; do
    fullpath="$SCRIPT_DIR/$script"

    if [ ! -f "$fullpath" ]; then
        echo "Skipping missing script: $script"
        continue
    fi

    if [ -z "$prev_jobid" ]; then
        jobid=$(sbatch "$fullpath" | awk '{print $4}')
    else
        jobid=$(sbatch --dependency=afterok:$prev_jobid "$fullpath" | awk '{print $4}')
    fi

    echo "Submitted $script as job $jobid"
    prev_jobid=$jobid

    sleep 1
done
