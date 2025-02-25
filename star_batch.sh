## Bach script for alignment using STAR
## Author: Simranjeet Kaur

#!/bin/bash

#PBS -W group_list=c******** -A c********
#PBS -N STAR_Alignment
#PBS -o /home/projects/c********/logs/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e /home/projects/c********/logs/$PBS_JOBNAME.$PBS_JOBID.err
#PBS -l nodes=2:ppn=30  # Adjust ppn based on available cores
#PBS -l walltime=999:00:00

# Load required modules
module load gcc/7.2.0
module load star/2.7.9a  # Load STAR module
module load parallel/20220422 # Load GNU Parallel module

# Define FASTQ directories
FASTQ_DIR_1="/home/projects/c******/fastq_clean"
FASTQ_DIR_2="/home/projects/c******/clean"

# Define STAR genome index and output directory
STAR_INDEX="/home/projects/c******/Reference/GRCh38/GRCh38_107/STAR_index"
OUTPUT_DIR="/home/projects/c******/star_alignment"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Export variables to make them accessible inside functions
export FASTQ_DIR_1 FASTQ_DIR_2 OUTPUT_DIR STAR_INDEX


# List of sample prefixes

SAMPLES=("S1-9-CTRL" "S5-10-CTRL" "S9-11-CTRL" "S13-12-CTRL" "S3-9-MIX" "S7-10-MIX" "S11-11-MIX" "S15-12-MIX" "HI-1A-19_1.CTRL" "HI-1A-19_2.IL+IFN" "HI-1A-19_3.IFNa" "HI-05C-19_1.CTRL" "HI-05C-19_2.IL+IFN" "HI-05C-19_3.IFNa" "HI-06C-19_1.CTRL" "HI-06C-19_2.IL+IFN" "HI-06C-19_3.IFNa" "HI-06D-18_1.CTRL" "HI-06D-18_2.IL+IFN" "HI-06D-18_3.IFNa" "HI-07B-19_1.CTRL" "HI-07B-19_2.IL+IFN" "HI-07B-19_3.IFNa")


# Function to find FASTQ files in both directories
##NOTE: DIR_2 has subfolders so do the following {FASTQ_DIR_2}/${SAMPLE}/${SAMPLE} to find the fastq files

find_fastq() {
    SAMPLE=$1
    if [ -f "${FASTQ_DIR_1}/${SAMPLE}_1.fq.gz" ]; then
        echo "${FASTQ_DIR_1}/${SAMPLE}_1.fq.gz,${FASTQ_DIR_1}/${SAMPLE}_2.fq.gz"
    elif [ -f "${FASTQ_DIR_2}/${SAMPLE}/${SAMPLE}_1.fq.gz" ]; then
        echo "${FASTQ_DIR_2}/${SAMPLE}/${SAMPLE}_1.fq.gz,${FASTQ_DIR_2}/${SAMPLE}/${SAMPLE}_2.fq.gz"
    else
        echo "Error: FASTQ files for sample $SAMPLE not found in either directory" >&2
        exit 1
    fi
} 

## If there are no subfolders do the following

#find_fastq() {
#    SAMPLE=$1
#    if [ -f "${FASTQ_DIR_1}/${SAMPLE}_1.fq.gz" ]; then
#        echo "${FASTQ_DIR_1}/${SAMPLE}_1.fq.gz,${FASTQ_DIR_1}/${SAMPLE}_2.fq.gz"
#    elif [ -f "${FASTQ_DIR_2}/${SAMPLE}_1.fq.gz" ]; then
#        echo "${FASTQ_DIR_2}/${SAMPLE}_1.fq.gz,${FASTQ_DIR_2}/${SAMPLE}_2.fq.gz"
#    else
#        echo "Error: FASTQ files for sample $SAMPLE not found in either directory" >&2
#        exit 1
#    fi
#}

# Function to run STAR for a single sample
run_star() {
    SAMPLE=$1
     echo "Running STAR alignment for sample: $SAMPLE"

    # Get the FASTQ file paths
    FASTQ_PAIR=$(find_fastq "$SAMPLE")
    FASTQ_1=$(echo "$FASTQ_PAIR" | cut -d',' -f1)
    FASTQ_2=$(echo "$FASTQ_PAIR" | cut -d',' -f2)

    # Define the sample-specific output directory
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
    
   # Create sample directory explicitly
    mkdir -p "$SAMPLE_OUTPUT_DIR" || { echo "Failed to create directory: $SAMPLE_OUTPUT_DIR"; exit 1; }

    # Run STAR
    STAR --runThreadN 4 \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$FASTQ_1" "$FASTQ_2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${SAMPLE_OUTPUT_DIR}/${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate

    echo "STAR alignment completed for sample: $SAMPLE"
}

export -f run_star  # Export the function for parallel execution
export -f find_fastq  # Export the helper function

# Run STAR in parallel for all samples
echo "Starting STAR alignment in parallel..."
parallel -j 7 run_star ::: "${SAMPLES[@]}"  # Adjust -j to the number of available cores

echo "All STAR alignments completed!"

