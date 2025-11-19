#!/bin/bash

#PBS -W group_list=cu_10031 -A cu_10031
#PBS -N Tophat_batch_run_p13
#PBS -o $PBS_JOBNAME.$PBS_JOBID.out
#PBS -e $PBS_JOBNAME.$PBS_JOBID.err
#PBS -l nodes=1:ppn=20
#PBS -l walltime=99:00:00

module load tophat/2.1.1
module load samtools/1.9
module load bowtie2/2.3.2

# Define the input and output directories
REFERENCE_DIR="/Reference/GRCh38/GRCh38_107"
REFERENCE_GENOME="${REFERENCE_DIR}/GRCh38.107_reordered"
GTF_FILE="${REFERENCE_DIR}/GRCh38_107_reordered.gtf"
FASTQ_DIR="clean/"
OUTPUT_DIR="HI_tophat_alignment/"


SAMPLES=("01-21-ICTRL" "01-21-II2CYT" "01-21-III3CYT" "02-21-ICTRL" "02-21-II2CYT" "02-21-III3CYT" "E-07-22_3CYT" "E-07-22_CTRL" "E-07-22_2CYT" "E-14-22_3CYT" "E-14-22_2CYT" "E-14-22_3CYT")


# Loop through each sample and run Tophat
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"

    # Construct the input FASTQ file paths
    FASTQ_FILES="${FASTQ_DIR}/${SAMPLE}_1.fq.gz,${FASTQ_DIR}/${SAMPLE}_2.fq.gz"

    # Construct the output directory for the current sample
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"

    # Run Tophat
    tophat -G "${GTF_FILE}" -p8 --library-type=fr-firststrand -o "${SAMPLE_OUTPUT_DIR}" \
           "${REFERENCE_GENOME}" "${FASTQ_FILES}"
done

echo "Batch processing completed!"

