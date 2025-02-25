## Bash script for creating read counts from sorted bam files using htseqcount using parallel processing
## Author: Simranjeet Kaur

#!/bin/bash

#PBS -W group_list=c******** -A c********/
#PBS -N HTSeq_Count_Parallel
#PBS -o /home/projects/c********/logs/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e /home/projects/c********/logs/$PBS_JOBNAME.$PBS_JOBID.err
#PBS -l nodes=1:ppn=16  # Adjust ppn based on available cores
#PBS -l walltime=24:00:00

# Load required modules
module load anaconda2/4.0.0
module load parallel/20220422  # Load GNU Parallel (or install if not available)

# Define directories
#INPUT_DIR="/home/projects/c********/star_alignment"
OUTPUT_DIR="/home/projects/c********/htseq_counts_p13_HI"
GTF_FILE="/home/projects/c********/Reference/GRCh38/GRCh38_107/GRCh38.107.reordered.gtf"


# Export variables to make them accessible inside functions
#export INPUT_DIR OUTPUT_DIR GTF_FILE
export OUTPUT_DIR GTF_FILE

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List of sample directories (add sample names here)
SAMPLES=("S1-9-CTRL" "S5-10-CTRL" "S9-11-CTRL" "S13-12-CTRL" "S3-9-MIX" "S7-10-MIX" "S11-11-MIX" "S15-12-MIX" "HI-1A-19_1.CTRL" "HI-1A-19_2.IL+IFN" "HI-1A-19_3.IFNa" "HI-05C-19_1.CTRL" "HI-05C-19_2.IL+IFN" "HI-05C-19_3.IFNa" "HI-06C-19_1.CTRL" "HI-06C-19_2.IL+IFN" "HI-06C-19_3.IFNa" "HI-06D-18_1.CTRL" "HI-06D-18_2.IL+IFN" "HI-06D-18_3.IFNa" "HI-07B-19_1.CTRL" "HI-07B-19_2.IL+IFN" "HI-07B-19_3.IFNa")


# Function to process a single sample
process_sample() {
    SAMPLE=$1
    echo "Processing sample: $SAMPLE"

    # Define input and output file paths
    SORTED_BAM_FILE="${OUTPUT_DIR}/sortedbyname_${SAMPLE}.bam"
    HTSEQ_OUTPUT_FILE="${OUTPUT_DIR}/geneIDCounts.${SAMPLE}_reverse.txt"

    # Run htseq-count
    echo "Running HTSeq-count for sample: $SAMPLE"
    htseq-count -f bam -s reverse -r name "$SORTED_BAM_FILE" "$GTF_FILE" > "$HTSEQ_OUTPUT_FILE"

    echo "Completed processing for sample: $SAMPLE"
}

export -f process_sample  # Export the function for GNU Parallel

# Run the function in parallel for all samples
parallel -j 4 process_sample ::: "${SAMPLES[@]}"  # Adjust -j to match the number of available cores

echo "All samples processed in parallel!"

