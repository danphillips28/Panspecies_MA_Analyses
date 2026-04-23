#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"

# List the scripts explicitly, in the order you want them submitted
scripts=(
	 	#"15_Downloading_C_elegans_TimeSeries_Expression2.sh"
	 	#"15_Downloading_Dmelanogaster_Epigenetics_Covariates_Data.sh"
	 	#"15_Downloading_D_melanogaster_TimeSeries_Expression.sh"
	 	#"15_Downloading_D_melanogaster_XRseq_Replicates_and_Merging.sh"
	 	#"15_Downloading_Hsapiens_Epigenetics_Covariates_Data.sh"
	 	#"15_Downloading_Ppacificus_Epigenetics_Covariates_Data.sh"
	 	#"15_Downloading_P_pacificus_TimeSeries_Expression_Submitter.sh")
	 	# Associating Mutation Rate and TimeSeries Expression
         	#"18_Associating_TimeSeries_Expression_with_Mutation_Rate_Submitter.sh")
	 	# Associating Mutation Rate with Expression
         	"16_Associating_Expression_with_Mutation_Rate_Submitter.sh"
	 	# Associating Mutation Rate and Epigenetics
	 	#"16_Associating_Epigenetics_and_Expression_etc_with_Mutation_Rate_Submitter.sh"
	 # Associating Change in Expression with Change in Mutation Rate Across C. elegans and P. pacificus
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
