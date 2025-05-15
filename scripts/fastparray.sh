#!/bin/bash
#SBATCH
#SBATCH --array=1-2%10
#SBATCH --partition=iob_p
#SBATCH --job-name=fastp
#SBATCH --job-name=fastp
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --mem=40GB
#SBATCH --output=./slurmoutputs/fastp/L001/fastp.%a.out
#SBATCH --error=./slurmoutputs/fastp/L001/fastp.%a.out
#SBATCH --mail-type=END

# Script to trim primer sequences from the amplicons using cutadapt
# Assumes we are in the working directory (./)
# Will make the necessary sub-directories if they are not present:
# ./trimmed${dataset}${run}/ # to hold trimmed sequence files
# ./slurmoutputs/cutadapt/ # to hold output overviews
# Assumes raw sequence data is held in a folder with the following structure:
# ./raw${run}

# Load fastp
module purge
module load fastp/0.23.2-GCC-11.3.0

# Make a directory to put the trimmed samples in
mkdir -p ./trimmed/ # Base directory for trimmed files
mkdir -p ./trimmed/${dataset}
mkdir -p ./fastpreports/
mkdir -p fastpreports/${dataset}
mkdir -p fastpreports/${dataset}/html
mkdir -p fastpreports/${dataset}/json

# Get input file names (paired forward/reverse reads)
FWD=`ls ./${inputdirectory}/ | grep "_R1_" | head -n $SLURM_ARRAY_TASK_ID | tail -1` # forward read
REV=`echo $FWD | sed 's/_R1_/_R2_/'` # paired reverse read

# Extract sample name from the full input filename:
#SAMPLE=$(echo ${FWD} | sed "s/_R1_\001\.fastq.gz*//")
SAMPLE=$(basename "${FWD}" | sed -E "s/.*_(S[0-9]+_L[0-9]+)_R1_001\.fastq\.gz/\1/")

echo "Beginning fastp processing of ${SAMPLE} ..."

fastp \
 -i ${inputdirectory/}/$FWD \
 -I ${inputdirectory}/$REV \
 -o ./trimmed/${dataset}/$FWD \
 -O ./trimmed/${dataset}/$REV \
 --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
 -f 30 \
 -F 30 \
 -t 5 \
 -T 5 \
 -g \
 -q 20 \
 -u 40 \
 -y \
 -Y 30 \
 -l 100 \
 -c \
 -p \
 -w 16 \
 -h "${SAMPLE}_fastp.html" \
 -j "${SAMPLE}_fastp.json"

echo "fastp processing of ${SAMPLE} complete..."

echo "fastp end..."

# Move all reports to the 'fastpreports' directory

echo "Moving fastp reports to final destination..."

mv *${SAMPLE}_fastp.html ./fastpreports/${dataset}/html/
mv *${SAMPLE}_fastp.json ./fastpreports/${dataset}/json/

echo "fastp quality filtering. Check the reports to verify outputs."
