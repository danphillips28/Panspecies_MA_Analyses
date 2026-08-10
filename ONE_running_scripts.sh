#!/bin/bash

#SBATCH --mail-user=daniel.phillips@lmh.ox.ac.uk
#SBATCH --mail-type=ALL

# Optional cleanup step
rm -rf /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data/
rm -rf /home/ocdm0351/DPhil/logs/
rm -rf /home/ocdm0351/DPhil/R_Data/
cp /home/ocdm0351/DPhil/scripts/misc/plotting.themes.R /home/ocdm0351/DPhil/R_Data/plotting.themes.R

SCRIPT_DIR="/home/ocdm0351/DPhil/scripts"

# List the scripts explicitly, in the order you want them submitted
scripts=(
"0_running_scripts.sh"
"1_Cleaning_Variant_Files.Rmd"
"1_Cleaning_Variant_Files_Submitter.sh"
"1_Removing_Bad_Samples_from_Behringer.sh"
"2_grooming_VCFs.sh"
"3_renaming_Ncrassa_contigs.sh"
"4_liftingOver.sh"
"5_building_snpEff_databases.sh"
"6_snpEffSift.sh"
"7_Downloading_Celegans_Essential_Genes.sh"
"7_Downloading_Celegans_HouseKeeping_Genes.sh"
"7_Downloading_Dmelanogaster_Essential_Genes.sh"
"8_Cleaning_Annotated_Variant_Files_Submitter.sh"
"9_Summarising_VCFs_Submitter.sh"
"10_EnrichmentAroundFeatures_Submitter.sh" # Need to run this because it produces some genomic resource objects needed in other scripts
) # Need to create sister species manually before running next section of pipeline # Currently automated at end of this script

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

cp /home/ocdm0351/DPhil/scripts/misc/Sister_Species_Table.txt /home/ocdm0351/DPhil/R_Data/Sister_Species_Table.txt
