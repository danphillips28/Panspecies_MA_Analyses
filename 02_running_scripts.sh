#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"

# List the scripts explicitly, in the order you want them submitted
scripts=(
	 # Associating Mutation Rate and Divergence
	 "10_Creating_Sister_Species_Table.sh"
	 "11_Renaming_References_for_OrthologR.sh"
	 "12_Calculating_dNdS_Submitter.sh"
	 "13_Associating_dNdS_with_Mutation_Rate_Submitter.sh"
	 "13_5_MutationRate_Functional_Enrichment_Submitter.sh"
	 # Associating Mutation Rate and Expression and Epigenetics
	 	# Downloading Expression and Covariates
	 	"14_Downloading_H_sapiens_Expression.sh"
	 	"15_Downloading_Celegans_Epigenetics_Covariates_Data.sh"
	 	"15_Downloading_C_elegans_H3K9me2_Replicates_and_Merging.sh"
	 	"15_Downloading_C_elegans_H3K9me3_Replicates_and_Merging.sh"
	 	"15_Downloading_C_elegans_HPL2_Replicates_and_Merging.sh")
		# Need to stop now and do some manual stuff at the start of 15_Downloading_C_elegans_TimeSeries_Expression2.sh

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
