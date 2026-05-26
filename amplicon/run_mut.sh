#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=matchbox

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=48:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

# SBATCH --array=1-2

# Set the file for output (stdout)
# SBATCH --output=stdout

# Set the file for error log (stderr)
# SBATCH --error=stderr

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# SBATCH --partition=genomics
# SBATCH --qos=genomics

# SBATCH --dependency=afterok:10567760
# SBATCH --reservation=highmem

# Command to run a serial job

# Create tmp directory for intermediate files
mkdir -p tmp

for file in P*/P*_r1.fq; do

dir=$(dirname $file)
base=$(basename $file _r1.fq)
r2=$(echo $file | sed 's/r1/r2/')

# first round exact match
# unmatched reads are saved in tmp/no_r1.fq and tmp/no_r2.fq
~/yan.a/software/matchbox/target/release/matchbox -s search_multi_pe.mb --error 0 --threads 16 $file --paired-with $r2 
mv mate_hits.tsv tmp/${base}_mate_hits.tsv

# second round with error 0.1 for srsf2
~/yan.a/software/matchbox/target/release/matchbox -s srsf2_10bp.mb --error 0.1 --threads 16 tmp/no_r1.fq
if [ ! -f second_hits.tsv ]; then
    touch second_hits.tsv
fi
mv second_hits.tsv tmp/${base}_second_R1_hits.tsv

~/yan.a/software/matchbox/target/release/matchbox -s srsf2_10bp.mb --error 0.1 --threads 16 tmp/no_r2.fq
if [ ! -f second_hits.tsv ]; then
    touch second_hits.tsv
fi
mv second_hits.tsv tmp/${base}_second_R2_hits.tsv

# third round with error 0.2 for srsf2
~/yan.a/software/matchbox/target/release/matchbox -s srsf2_10bp.mb --error 0.2 --threads 16 tmp/no_r1.fq
if [ ! -f second_hits.tsv ]; then
    touch second_hits.tsv
fi
mv second_hits.tsv tmp/${base}_third_R1_hits.tsv

~/yan.a/software/matchbox/target/release/matchbox -s srsf2_10bp.mb --error 0.2 --threads 16 tmp/no_r2.fq
if [ ! -f second_hits.tsv ]; then
    touch second_hits.tsv
fi
mv second_hits.tsv tmp/${base}_third_R2_hits.tsv

# Merge all TSV files and modify second column based on filename
echo "Merging and processing TSV files..."

# Process and merge all TSV files
for tsv_file in tmp/*_hits.tsv; do
    if [ -s "$tsv_file" ]; then  # Only process non-empty files
        filename=$(basename "$tsv_file")
        
        # Determine the new second column value based on filename
        if [[ $filename == *"_mate_hits.tsv" ]]; then
            # For mate_hits, keep original second column (R1 or R2rc)
            cat "$tsv_file"
        elif [[ $filename == *"_second_R1_hits.tsv" ]]; then
            # Replace second column with R1_2
            awk -F'\t' 'BEGIN{OFS="\t"} {if(NF>=2) $2="R1_2"; print}' "$tsv_file"
        elif [[ $filename == *"_second_R2_hits.tsv" ]]; then
            # Replace second column with R2rc_2
            awk -F'\t' 'BEGIN{OFS="\t"} {if(NF>=2) $2="R2rc_2"; print}' "$tsv_file"
        elif [[ $filename == *"_third_R1_hits.tsv" ]]; then
            # Replace second column with R1_3
            awk -F'\t' 'BEGIN{OFS="\t"} {if(NF>=2) $2="R1_3"; print}' "$tsv_file"
        elif [[ $filename == *"_third_R2_hits.tsv" ]]; then
            # Replace second column with R2rc_3
            awk -F'\t' 'BEGIN{OFS="\t"} {if(NF>=2) $2="R2rc_3"; print}' "$tsv_file"
        fi
    fi
done > ${dir}/${base}_merged_all_hits.tsv

rm tmp/*

done

