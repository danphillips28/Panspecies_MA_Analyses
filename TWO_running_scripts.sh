#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"

# List the scripts explicitly, in the order you want them submitted
scripts=(
#"11_Creating_Sister_Species_Table.sh"
#"12_Renaming_References_for_OrthologR.sh"
#"13_Calculating_dNdS_Submitter.sh"
#"14_Associating_dNdS_with_Mutation_Rate_Submitter.sh"
	# Downloading Expression
#"15_Downloading_C_elegans_Expression_Submitter.sh"
#"15_Downloading_D_melanogaster_Expression.sh"
#"15_Downloading_P_pacificus_Expression_Submitter.sh"
	# Functional Enrichment
#"16_MutationRate_Functional_Enrichment_Submitter.sh"
	# Downloading C elegans Covariates
#"17_Downloading_C_elegans_TimeSeries_Expression2"
#"17_Downloading_Celegans_Epigenetics_Covariates_Data.sh"
#"17_Downloading_C_elegans_H3K9me2_Replicates_and_Merging.sh"
#"17_Downloading_C_elegans_H3K9me3_Replicates_and_Merging.sh"
#"17_Downloading_C_elegans_HPL2_Replicates_and_Merging.sh"
	# Downloading D melanogaster Covariates
#"17_Downloading_Dmelanogaster_Epigenetics_Covariates_Data.sh"
#"17_Downloading_D_melanogaster_TimeSeries_Expression.sh"
#"17_Downloading_D_melanogaster_XRseq_Replicates_and_Merging.sh"
	# Downloading P pacificus Covariates
"17_Downloading_Ppacificus_Epigenetics_Covariates_Data.sh"
"17_Downloading_P_pacificus_TimeSeries_Expression_Submitter.sh"
	# Downloading H sapiens Data
#"17_Downloading_Hsapiens_Epigenetics_Covariates_Data.sh"
#"17_Downloading_H_sapiens_Expression.sh"
)
# Need to stop now and do some manual stuff for Cel timeseries - it's the liftOver. Only needed when running on new location, i have files saved.

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
