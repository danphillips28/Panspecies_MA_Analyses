#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

# Optional cleanup step
rm -rf /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data/
rm -rf /home/ocdm0351/DPhil/logs/

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"

# List the scripts explicitly, in the order you want them submitted
scripts=(
	 # Some Formattng and Grooming
	 "1_Cleaning_Variant_Files_Submitter.sh" # submits 1_Cleaning_Variant_Files.Rmd
	 "1_Removing_Bad_Samples_from_Behringer.sh"
	 "2_grooming_VCFs.sh"
	 "3_renaming_Ncrassa_contigs.sh"
	 # liftOver of Variants if Needed
	 "4_liftingOver.sh"
	 # Variant Annotation
	 "5_building_snpEff_databases.sh"
	 "6_snpEffSift.sh"
	 "6_5_Downloading_Celegans_Essential_Genes.sh"  
	 "6_5_Downloading_Celegans_HouseKeeping_Genes.sh"  
	 "6_5_Downloading_Dmelanogaster_Essential_Genes.sh"
	 "7_Cleaning_Annotated_Variant_Files_Submitter.sh"
	 # Summary Statistics of Variants
	 "8_Summarising_VCFs_Submitter.sh"
	 # Testing Enrichment of Variants Around Gene Bodies
	 "9_EnrichmentAroundFeatures_Submitter.sh") # Need to manually create Sister_Species file now before running script #10

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
