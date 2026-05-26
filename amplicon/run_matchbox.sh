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


# Get unique P numbers from barcode.tsv
plateid=$(awk -F'\t' 'NR>1 {print $1}' output/barcode_renamed.tsv | cut -d'_' -f1 | sort | uniq)

# Split barcode.tsv into subsets based on P number
for p in $plateid; do

    mkdir -p ${p}

    cd ${p}
    # Get header
    head -n 1 ../output/barcode_renamed.tsv > ${p}_barcode.tsv
    # Get rows for this P number
    grep "^${p}" ../output/barcode_renamed.tsv >> ${p}_barcode.tsv

    # Find all fastq files for this plate
    pid=$(echo ${p} | sed 's/^P//g')
    # fastq_files=$(ls ../*-${p}-*bbmerged.fastq.gz  2>/dev/null)
    fastq_files1=$(ls ../../fastq/NGU17315r/Plate-${pid}*_R1_001.fastq.gz  2>/dev/null)
    fastq_files2=$(ls ../../fastq/NGU17315r/Plate-${pid}*_R2_001.fastq.gz  2>/dev/null)

    if [ -n "$fastq_files1" ]; then
        echo "Processing plate ${p}..."
        # Create matchbox script for this plate
        cat > ${p}_demux.mb << EOF

bcs = tsv('${p}_barcode.tsv')   

if read.r1.concat(-read.r2) is [ start:|1| bc.forward _ -(bc.reverse) end:|1| ] for bc in bcs => {
    count!(bc.name)
    read.r1.out!('{bc.name}_r1.fq')
    read.r2.out!('{bc.name}_r2.fq')
}
[_] => {
    count!('unmatched') 
    read.r1.out!('unmatched_r1.fq') 
    read.r2.out!('unmatched_r2.fq')   
}
EOF
        # Run matchbox for this plate
        echo "Running matchbox for plate ${p} using ${fastq_files1} and ${fastq_files2} and ${p}_demux.mb"
        ~/yan.a/software/matchbox/target/release/matchbox -s ${p}_demux.mb --error 0.15 --threads 16 $fastq_files1 --paired-with $fastq_files2
    else
        echo "No fastq files found for plate ${p}"
    fi

    cd ..

done

    
    
