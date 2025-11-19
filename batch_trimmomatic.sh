#!/bin/bash

# Define paths to input files, adapter sequence, and output directory
INPUT_DIR="/path/to/your/input_data"
OUTPUT_DIR="/path/to/your/output_data"
ADAPTER="CTGTCTCTTATACACATCT"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# List all fastq files in the input directory (paired-end files assumed)
for f in $INPUT_DIR/*_R1_*.fastq.gz; do
    # Get the base name for the files
    base=$(basename $f "_R1_001.fastq.gz")

    # Define the corresponding R2 file (reverse read)
    R1="$INPUT_DIR/${base}_R1_001.fastq.gz"
    R2="$INPUT_DIR/${base}_R2_001.fastq.gz"

    # Define the output file names
    OUTPUT_R1="$OUTPUT_DIR/${base}_R1_trimmed.fastq.gz"
    OUTPUT_R2="$OUTPUT_DIR/${base}_R2_trimmed.fastq.gz"
    OUTPUT_UNPAIRED_R1="$OUTPUT_DIR/${base}_R1_unpaired.fastq.gz"
    OUTPUT_UNPAIRED_R2="$OUTPUT_DIR/${base}_R2_unpaired.fastq.gz"

    # Run Trimmomatic in parallel
    echo "Trimming $base"
    trimmomatic PE -threads 4 \
        $R1 $R2 \
        $OUTPUT_R1 $OUTPUT_UNPAIRED_R1 \
        $OUTPUT_R2 $OUTPUT_UNPAIRED_R2 \
        ILLUMINACLIP:$ADAPTER:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36 &
done

echo "Trimming completed for all samples!"
