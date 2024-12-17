### ****************************** README ********************************** ###

################################################################################
#################### Introduction and Important Notes ##########################
################################################################################

# This script was compiled by Philip Brailey-Jones for the Bennetzen Lab
# October 2023
# Last Update: 12th August 2024

# Script assumes that all sbatch scripts are submitted in the working directory (./)

################################################################################
############################ Required Programs #################################
################################################################################

#**************** Shotgun WGS Pipeline ****************#

#Programs available on GACRC Sapelo2 updated as of October 2023

#fastp
fastp/0.23.2-GCC-11.3.0

#************* 16S gRNA Amplicon Pipeline *************#

################################################################################
############################ Required R Packages ###############################
################################################################################

################################################################################
############################# Directory Structure ##############################
################################################################################

# All paths are relative to the working directory (.)
# The majority of these directories will be made as part of the pipeline so you do not necessarily need to make them yourself here

# Directory for all scripts
mkdir -p ./scripts/

# Directory for all output and error files
# Some steps where parralelization is used will create additional subdirectories themselves to not fill up this directory
mkdir ./slurmoutputs

# Directory of reference databases used
mkdir -p ./databases
mkdir -p ./databases/sorghumv3/

# Directory for raw sequence input data:
# You WILL need to make and populate this directory with your raw files
mkdir -p ./raw/

# Directory for primer-trimmed sequence data from cutadapt
mkdir -p ./trimmed/

# Directory for fastQC output files:
mkdir -p ./fastqcoutput/










# Subdirectory for raw WGS sequence data:
mkdir ./raw_WGS/

# Subdirectory for sequences trimmed using fastp
mkdir ./trimmed_WGS/

# Subdirectory for fastp reports
mkdir ./fastpreports_WGS/

# Subdirectory for raw 16S gDNA amplicon sequence data:
mkdir ./raw_16S/

# Subdirectory for cutada[t trimmed 16S gDNA amplicon sequence data:
mkdir ./trimmed_16S/

#Subdirectory to download all databases to
mkdir ./databases

	# Sub-subdirectory containing
	mkdir ./databases/sorghumv3

	# Sub-subdirectory containing the WoL2 database
	mkdir ./databases/WoL2/

mkdir bowtie2outputs_WGS

mkdir bowtie2outputs_WGS/outputs

mkdir bowtie2outputs_WGS/hostalignedonly

mkdir bowtie2outputs_WGS/unalignedonly

mkdir bowtie2outputs_WGS/WoL2outputs

mkdir bowtie2outputs_WGS/WoL2alignedonly

# Final directories

mkdir ./databases/greengenes2
mkdir ./greengenes2WGS/


wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.nwk.qza
mv 2022.10.taxonomy.asv.nwk.qza databases/greengenes2

########################################################
############# Required files and locations #############
########################################################

# All paths are relative to the working directory (.)


########################################################
########################################################
########################################################
##################### WGSS Pipeline ####################
########################################################
########################################################
########################################################

# Download sorghum reference genome
# We will download the v.3 genome from NCBI. I looked into downloading the v.5 genome from phytozome but they just pointed to the v.3 genome for some reason (07/12/2023)

# I couldn't find the ftp link for chloroplast or mitochondria so I just copy/pasted them into a fasta titled 'host_plastids.fasta'- Thanks NCBI for trying to make everything fancier and hiding / removing hyperlinks!

# Submit batch script to download sorghum reference database
sbatch scripts/sorghumdl.sh

#------------------------ sorghumdl.sh START------------------------------------
#!/bin/bash
#SBATCH --job-name=sorghumdl
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=./slurmOutputs/sorghumdl.out
#SBATCH --error=./slurmOutputs/sorghumdl.out
#SBATCH --mail-user=pab14613@uga.edu
#SBATCH --mail-type=START,END

# Script to download the sorghum reference database
# Assumes we are in the working directory (./)
# Will make the necessary sub-directories if they are not present:
# ./databases/ # to hold all databases
# ./databases/sorghumv2 # to hold the sorghum genome

# Make required databases if needed

mkdir -p ./databases

mkdir -p ./databases/sorghumv3

# URL of the gzipped file

file_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz"

# Output file name (you can change this if needed)

output_file="host_genome.fasta"


# Download the gzipped file

echo "Downloading genome sequence..."

curl -O $file_url


# Unzip the downloaded file

echo "Downloading complete, unzipping..."

gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz > $output_file

# Move the sorghum genome reference to the intended final location

echo "Moving file to final destination..."

mv ./$output_file ./databases/sorghumv3/


# Optional: Remove the gzipped file if you don't need it anymore

rm GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz

echo "Genome sequence saved to $output_file"

#-----------------------SLURM SCRIPT sorghumdl.sh END---------------------------

#--------------------------- switchgrassdl.sh ----------------------------------
#!/bin/bash
#SBATCH --job-name=sorghumdl
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=./slurmOutputs/switchgrassdl.out
#SBATCH --error=./slurmOutputs/switchgrassdl.out
#SBATCH --mail-user=pab14613@uga.edu
#SBATCH --mail-type=START,END

# Script to download the sorghum reference database
# Assumes we are in the working directory (./)
# Will make the necessary sub-directories if they are not present:
# ./databases/ # to hold all databases
# ./databases/sorghumv2 # to hold the sorghum genome

# Make required databases if needed

mkdir -p ./databases

mkdir -p ./databases/switchgrassv5

# URL of the gzipped file

file_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/808/335/GCF_016808335.1_P.virgatum_v5/GCF_016808335.1_P.virgatum_v5_genomic.fna.gz"


# Output file name (you can change this if needed)

output_file="host_genome.fasta"


# Download the gzipped file

echo "Downloading genome sequence..."

curl -O $file_url


# Unzip the downloaded file

echo "Downloading complete, unzipping..."

gunzip -c GCF_016808335.1_P.virgatum_v5_genomic.fna.gz > $output_file

# Move the sorghum genome reference to the intended final location

echo "Moving file to final destination..."

mv ./$output_file ./databases/switchgrassv5/


# Optional: Remove the gzipped file if you don't need it anymore

rm GCF_016808335.1_P.virgatum_v5_genomic.fna.gz

echo "Genome sequence saved to $output_file"

#-----------------------SLURM SCRIPT sorghumdl.sh END---------------------------

#******************************************************************************#
#*********************** Sequence Quality Control *****************************#
#******************************************************************************#

# Fastp will generate quality reports both before and after quality filtering the sequences. In this iteration fastp automatically detects and removes adaptersls.
# We also filter out any sequences with a phred quality >=Q15; >40$ bases unqualified; length < 50bp.

# Run fastp.sh script
sbatch ./scripts/fastp.sh

#------------------------------fastp.sh START-----------------------------------

#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=./slurmoutputs/fastp.out
#SBATCH --error=./slurmoutputs/fastp.out
#SBATCH --mail-user=pab14613@uga.edu
#SBATCH --mail-type=END

# Script to perform quality control and trimming of raw illumina sequences
# Assumes we are in the working directory (./)
# Will make the necessary sub-directories if they are not present:
# ./trimmed/ # to hold all trimmed sequences
# ./fastpreports/ # To hold all fastp quality reports
# ./multiqcreport # To hold the aggregate multiqc report

echo "Making directories..."

# make necessary directories
mkdir -p ./trimmed
mkdir -p ./fastpreports
mkdir -p ./multiqcreport

echo "Loading modules..."

# Load fastp
module purge
module load fastp/0.23.2-GCC-11.3.0
module load MultiQC/1.14-foss-2022a

# Loop fastp quality filtering over both fwd and rev reads of each sample

echo "fastp start..."

cd ./raw/

for f1 in *_R1_001.fastq.gz

do

echo "Beginning fastp processing of ${f1} ..."

f2=${f1/R1/R2}

fastp \
 -i $f1 \
 -I $f2 \
 -o ../trimmed/$f1 \
 -O ../trimmed/$f2 \
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
 -h ${f1/%_R1*/}"_fastp_report.html" \
 -j ${f1/%_R1*/}"_fastp_report.json"

echo "fastp processing of ${f1} complete..."

done

echo "fastp end..."

# Move all reports to the 'fastpreports' directory

echo "Moving fastp reports to final destination..."

mv *report.html *report.json ../fastpreports/

# Create a multiQC report aggregating all samples
# SEP 24: This isn't aggregating the reports... needs troubleshooting.

#echo "Performing multiQC aggregation of fastp reports..."

#multiqc --force --outdir ../multiqcreport/ ../fastpreports/

echo "fastp quality filtering complete. Check the reports to verify outputs."

#-----------------------------=-fastp.sh END------------------------------------

# We can now review the fastp individual reports and the multiQC combined report

# Let's count how many sequences we have per sample before and after quality removal, we will attach these all together later

#raw input
sbatch --export=inputdirectory=./raw,outputfile=countsummary_rawseqs scripts/countreads.sh
#trimmed / quality filtered output
sbatch --export=inputdirectory=./trimmed,outputfile=countsummary_trimmedseqs scripts/countreads.sh

#----------------------------countreads.sh START--------------------------------

#!/bin/bash
#SBATCH --job-name=countreads
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=250G
#SBATCH --time=4:00:00
#SBATCH --output=./slurmOutputs/countreads.out
#SBATCH --error=./slurmOutputs/countreads.out
#SBATCH --mail-user=pab14613@uga.edu
#SBATCH --mail-type=END

# Script to count the number of sequences in a FASTQ formatted sequence file
# Assumes we are in the working directory (./)
# The USER provides the directory location for files that will be counted
#       - This script will count all files in a directory indiscriminately. If you don't want them counting don't include them in the directory!

# Set the directory where your 'fastq.gz' files are located

input_dir="${inputdirectory}"

# Define the output file for the summary table

output_file="${outputfile}.csv"

# Create a header for the output CSV file

echo "Sample Name,Read Count" > "$output_file"

# Iterate over the 'fastq.gz' files in the directory

for file in "$input_dir"/*_R1_*.fastq.gz; do

	# Extract the sample name from the file name

	sample_name=$(basename "$file" | sed 's/\.fastq\.gz//')


	# Use zcat to extract the contents of the gzipped file and count the lines (reads), subtract 1 to account for the header

	read_count=$(zcat "$file" | wc -l)

	# Calculate the actual read count

	read_count=$((read_count / 4))

	# We divide by 4 because one fastq read comprises of 4 lines: 1) identifier, 2) sequence, 3) blank ("starts with a '+', may contain the same info as line 1", 4) sequence quality scores)

	# Append the sample name and read count to the output CSV file

	echo "$sample_name,$read_count" >> "$output_file"

done

echo "Read count summary has been saved to $output_file"

#-----------------------------countreads.sh END---------------------------------

#******************************************************************************#
#***************** Host Sequence Removal / Subsetting *************************#
#******************************************************************************#

# We remove host reads by aligning the whole dataset against the reference genome of the host (in this case Sorghum bicolor), and subsetting the datasets to contain only reads that *do* (for a sorghum specific file) and *do not* (for everything else) align to the host.
# We ALSO first remove chloroplast and mitochondrial reads, also through alignment, to ensure that these are not possibly incorporated into our whole sorghum genome alignment
#    - This is unlikely... but it might also go a good way ytowards discriminating what of the unknown sequences are truly unknown versus unassigned 'plant origin' sequences

# We follow three overarching steps per the following tool: https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
# Step 1: bowtie2- mapping of (all) reads against host genome
# Step 2: samtools view- 'filter-flags' extraction of unmapped reads
# Step 3: samtools fastq- split paired-end reads into separated R1 and R2 fastq that can be used downstream

# Splitting reads into separated R1 and R2 files allows them to be put into subsequent steps effectively as paired forward and reverse reads

# Create a bowtie2 database using the reference sorghum genom
sbatch --export=refgenome=./databases/sorghumv3/host_genome.fasta,indexlocation=./databases/bowtie2/sorghumv3bowtie2db,dbname=host_DB scripts/bowtie2build.sh
# Create a bowtie2 database using the reference plastid genomes
sbatch --export=refgenome=databases/plastiddb/host_plastids.fasta,indexlocation=./databases/bowtie2/plastidbowtie2db,dbname=plastid_DB scripts/bowtie2build.sh

# switchgrass
sbatch --export=refgenome=./databases/switchgrassv5/host_genome.fasta,indexlocation=./databases/bowtie2/switchgrassv5bowtie2db,dbname=host_DB scripts/bowtie2build.sh


#-------------------SLURM SCRIPT bowtie2build.sh START -------------------------

#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=buildbowtie2db
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --output=./slurmoutputs/bowtie2build.out
#SBATCH --error=./slurmoutputs/bowtie2build.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END

# Script to create a bowtie2 database to map experimental sequences against
# Assumes we are in the working directory (./)
# Will make the necessary sub-directories if they are not present
#       - ./databases/bowtie2   # This holds all of the bowtie2 index folders
#       - We made the database directory in a previous step, we now make subdirectories for the bowtie2 specific database which we call in our ${index} command

# load bowtie2
# module versions available on the UGA GACRC Server updated OCT 2023
module purge
module load Bowtie2/2.4.5-GCC-11.3.0
module list

# Make bowtie2 database
mkdir -p ./databases/bowtie2
mkdir -p ${indexlocation}

# Build the bowtie2 database files
bowtie2-build ${refgenome}  $dbname --seed 42 -f --threads 16

# Extract database information
bowtie2-inspect -n $dbname

# Move constructed database files to where they need to be
mv *.bt2 ${indexlocation}

#---------------------- bowtie2build.sh END ------------------------------------

# Run bowtie2 mapping against the reference dataset in parallel

# IMPORTANT: Update 'nSamples' to specify the number of samples for parallelization to be split across.
# Count samples with the following
cd ./trimmed/ # Move into the directory containing the sequences
ls -1 | wc -l
cd ../ # Move back into the working directoryinput

# Define the number of samples
nSamples=30 # Quant experiment
nSamples=12 # Mock experiment

# First we carry out plastid alignment of the quality filtered sequences
sbatch --array=1-$nSamples%10 --export=referencedb=databases/bowtie2/plastidbowtie2db/plastid_DB,outdir=plastidalignment,inputdir=trimmed,sensitivity=--sensitive,mininsvalue=200,maxinsvalue=6000 scripts/bowtie2mapping.sh
# Following this initial filtering we'll carry out further alignment against the sorghum reference genome, but ONLY for the sequences that DID NOT align against the plastids
sbatch --array=1-$nSamples%10 --export=referencedb=databases/bowtie2/sorghumv3bowtie2db/host_DB,outdir=hostalignment,inputdir=alignments/plastidalignment/unmappedreads,sensitivity=--sensitive,mininsvalue=200,maxinsvalue=6000 scripts/bowtie2mapping.sh
sbatch --array=1-$nSamples%10 --export=referencedb=databases/bowtie2/WoL2/WoLr2,outdir=WoL2alignment,inputdir=alignments/hostalignment/unmappedreads,sensitivity=--very-sensitive,mininsvalue=200,maxinsvalue=6000 scripts/bowtie2mapping.sh


# Switchgrass project
sbatch --array=1-$nSamples%10 --export=referencedb=databases/bowtie2/switchgrassv5bowtie2db/host_DB,outdir=hostalignment,inputdir=trimmed,sensitivity=--sensitive,mininsvalue=200,maxinsvalue=6000 scripts/bowtie2mapping.sh

#----------------------- bowtie2mapping.sh START ----------------------------------

#!/bin/bash
#SBATCH
#SBATCH --array=1-2%10
#SBATCH --partition=iob_p
#SBATCH --job-name=bowtie2mapping
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=10
#SBATCH --time=4:00:00
#SBATCH --mem=100GB
#SBATCH --output=./slurmoutputs/bowtie2/mapping.%a.out
#SBATCH --error=./slurmoutputs/bowtie2/mapping.%a.out
#SBATCH --mail-type=END

# Script to map experimental sequences against a bowtie2-formatted database of experimental sequences
# Assumes we are in the working directory (./)
# Will make the necessary sub-directories if they are not present:

# Make required directories
mkdir -p ./alignments
mkdir -p ./alignments/$outdir/
mkdir -p ./alignments/$outdir/unmappedreads/
mkdir -p ./alignments/$outdir/mappedreads/

# load bowtie2 and samtools
# module versions available on the UGA GACRC Server updated OCT 2023
module purge
module load Bowtie2/2.4.5-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

echo "Preparing sample ID for alignment..."

# get input files (paired forward/reverse reads)
FWD=`ls ./${inputdir}/ | grep "_R1_" | head -n $SLURM_ARRAY_TASK_ID | tail -1` # get forward read
REV=`echo $FWD | sed 's/_R1_/_R2_/'` # get paired reverse read

echo "Sample IDs are..."
echo "$FWD"
echo "$REV"
echo "./${inputdir}/${FWD}"
echo "./${inputdir}/${REV}"

# extract sample name from the full input filename:
SAMPLE=$(echo ${FWD} | sed "s/_R1_\001\.fastq.gz*//")

# ********************* Step 1: Bowtie2 mapping ****************************** #

# map sequences against the host-sequence database, keep both aligned and unaligned reads (paired-ends reads)

echo "Beginning alignment of $SAMPLE against reference database..."

# This makes the following sense for layout readability but for some reason the inout script would prefer it all in one line...
#bowtie2 -p 12 -x ./${referencedb} \
#  $sensitivity \
#	-I ${mininsvalue} \ # A little smaller than our expected minimum fragment size
#	-X ${maxinsvalue} \ # Longer than our expected largest fragment size because there is a pretty decent 'tail' beyond our main large fragment ~1250bp and I want to maintain as many reads as possible after mapping
#	-1 ./trimmed/$FWD \
#	-2 ./trimmed/$REV \
#	-S ./alignments/$outdir/${SAMPLE}_mappedandunmapped.sam

# Actual input
bowtie2 -p 12 -x ./${referencedb} $sensitivity -I ${mininsvalue} -X ${maxinsvalue} -1 ./trimmed/$FWD -2 ./trimmed/$REV -S ./alignments/$outdir/${SAMPLE}_mappedandunmapped.sam

echo "Alignment complete, beginning post-alignment dataset filtering..."

# disable mixed mode

# convert .sam file to .bam file

samtools view -bS ./alignments/$outdir/${SAMPLE}_mappedandunmapped.sam > ./alignments/$outdir/${SAMPLE}_mappedandunmapped.bam

# ********************** Step 2: Samtools filtering ************************** #

# filter only unmapped reads (requires bam input)

echo "Filtering only unmapped reads..."

samtools view -b -f 12 -F 256 ./alignments/$outdir/${SAMPLE}_mappedandunmapped.bam > ./alignments/$outdir/${SAMPLE}_bothreadsunmapped.bam

echo "Filtering only mapped reads..."

# filter only mapped reads and output a bam file

samtools view -b -f 3 ./alignments/$outdir/${SAMPLE}_mappedandunmapped.bam > ./alignments/$outdir/${SAMPLE}_bothreadsmapped.bam

# filter only mapped reads and output a sam file

samtools view -h -f 3 ./alignments/$outdir/${SAMPLE}_mappedandunmapped.bam > ./alignments/$outdir/${SAMPLE}_bothreadsmapped.sam

# ***** Step 3: Split paired-end reads ******* #

echo "Splitting reads belonging to mapped or unmapped datasets back into forward and reverse for further processing..."

# unmapped reads only

samtools sort -n -m 5G -@ 2 ./alignments/$outdir/${SAMPLE}_bothreadsunmapped.bam -o ./alignments/$outdir/${SAMPLE}_bothreadsunmapped_sorted.bam

samtools fastq -@ 8 ./alignments/$outdir/${SAMPLE}_bothreadsunmapped_sorted.bam -1 ./alignments/$outdir/unmappedreads/${SAMPLE}_R1_001.fastq.gz -2 ./alignments/$outdir/unmappedreads/${SAMPLE}_R2_001.fastq.gz -0 /dev/null -s /dev/null -n

# mapped reads only

samtools sort -n -m 5G -@ 2 ./alignments/${outdir}/${SAMPLE}_bothreadsmapped.bam -o ./alignments/${outdir}/${SAMPLE}_bothreadsmapped_sorted.bam

samtools fastq -@ 8 ./alignments/${outdir}/${SAMPLE}_bothreadsmapped_sorted.bam -1 ./alignments/${outdir}/mappedreads/${SAMPLE}_R1_001.fastq.gz -2 ./alignments/${outdir}/mappedreads/${SAMPLE}_R2_001.fastq.gz -0 /dev/null -s /dev/null -n

echo "Alignment and filtering complete"

#-----------------------bowtie2mapping.sh END-----------------------------------

# Now we want to know how many reads were preserved at each of our plant mapping steps, so we will count the PLASTID MAPPED, SORGHUM MAPPED, and SORGHUM-NON-MAPPED READS
# Counting the reads that *did not* map to Sorghum will tell us how many input reads we had in our subsequent mappings
# Hopefully by the law of maths this should all add back up to original numbers...

# Plastid mapped
sbatch --export=inputdirectory=./alignments/plastidalignment/mappedreads,outputfile=countsummary_plastidmapped scripts/countreads.sh
# Plastid non-mapped
sbatch --export=inputdirectory=./alignments/plastidalignment/unmappedreads,outputfile=countsummary_plastidunmapped scripts/countreads.sh
# Sorghum mapped
sbatch --export=inputdirectory=./alignments/hostalignment/mappedreads,outputfile=countsummary_sorghummapped scripts/countreads.sh
# Sorghum non-mapped
sbatch --export=inputdirectory=./alignments/hostalignment/unmappedreads,outputfile=countsummary_sorghumunmapped scripts/countreads.sh

# Switchgrass project
# Sorghum mapped
sbatch --export=inputdirectory=./alignments/hostalignment/mappedreads,outputfile=countsummary_switchgrassmapped scripts/countreads.sh
# Sorghum non-mapped
sbatch --export=inputdirectory=./alignments/hostalignment/unmappedreads,outputfile=countsummary_switchgrassunmapped scripts/countreads.sh



#******************************************************#
#************** Host Read Normalization ***************#
#******************************************************#

# we now have a file of reads that mapped to the sorghum genome- we want to count how many reads of the sorghum genome are present in each file. We will use this value for normalization later.

# Organize samples in a table with read counts, and normalization factors.
# Normalization factor is determined by the correction factor required for each individual samples' read count to be equal to the lowest samples' read count.
# We normalize in this way because the same mass of plant sample was used for each original DNA extraction.

# Create a bowtie2 database using the reference plastid genomes



#******************************************************************************#
#************************** Taxonomy Assignment *******************************#
#******************************************************************************#

# So we have lots of fun cleaned-up sequences to work with now, if only we can figure out what they are...



# The authors of Greengenes2 have integrated WGS metagenomics classification with 16S rDNA classification through the use of a unified database that classifications can be drawn on.
# To use the Greengenes2 database on WGS data we must generate Operational Genomic Units (OGUs) through Woltka. We can then place these OGUs into the tree-based taxonomy employed by Greengenes2.

#*************************** Woltka OGU taxonomy ******************************#

# Download the current iteration of the Web of Life 2 database

sbatch --export=dbformat=bowtie2,outdir=WoL2,urlinput=http://ftp.microbio.me/pub/wol2/databases/bowtie2/ ./scripts/dbdownload.sh
sbatch --export=dbformat=kraken2,outdir=WoL2,urlinput=http://ftp.microbio.me/pub/wol2/databases/kraken2/ ./scripts/dbdownload.sh
sbatch --export=dbformat=bracken,outdir=WoL2,urlinput=http://ftp.microbio.me/pub/wol2/databases/bracken/ ./scripts/dbdownload.sh

------------BASH SCRIPT WoL2download.sh START--------------

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=WoL2download
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=300G
#SBATCH --output=./slurmoutputs/dbdownload.out
#SBATCH --error=./slurmoutputs/dbdownload.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END

# Make the directory where we will download the files
mkdir -p ./databases/${dbformat}
mkdir -p ./databases/${dbformat}/${outdir}/

#******************************* db download **********************************#

echo "Downloading the ${dbformat} formatted ${outdir} database"

# URL of the database file
URL="${urlinput}"

# Directory where you want to save the downloaded file
DOWNLOAD_DIR="./databases/${dbformat}/${outdir}/"

# Download the file
wget -P "$DOWNLOAD_DIR" -r --no-parent -nd --cut-dirs=3 -A '*' "$URL"

echo "Download complete."

-------------BASH SCRIPT WoL2download.sh END---------------

# Align samples against the WoL2 database using bowtie2
# Following alignment we once again use samtools to subset this time only the mapped reads which we can then use in Woltka to generate OGUs.

# We do not need to make a bowtie2 format database this time- We have already downloaded the WoL2 bowtie formatted database at './databases/WoL2/'

nSamples=30 # Update this to match your number of samples (1/2 number of files in the file directory if you have FWD and REV reads)
# We can call back to our alignment script for the WoL2 alignment
sbatch --array=1-$nSamples%10 --export=referencedb=./databases/bowtie2/WoL2/WoLr2,outdir=WoL2alignment,inputdir=alignments/hostalignment/unmappedreads,sensitivity=--very-sensitive,mininsvalue=200,maxinsvalue=6000 scripts/bowtie2mapping.sh

# Count the mapped and unmapped reads
# WoL2 mapped
sbatch --export=inputdirectory=./alignments/WoL2alignment/mappedreads,outputfile=countsummary_wol2mmapped scripts/countreads.sh
# WoL2 unmapped
sbatch --export=inputdirectory=./alignments/WoL2alignment/unmappedreads,outputfile=countsummary_wol2unmapped scripts/countreads.sh

# This script generates a folder full of alignments against the WoL database
# You can now use these alignments as inputs to the OGU generation through Woltka

sbatch --export=indir=WoL2alignment ./scripts/woltkaOGU.sh

--------------SBATCH SCRIPT WoltkaOGU.sh START-------------

#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=woltkaogu
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=500G
#SBATCH --output=./slurmoutputs/woltkaogu.out
#SBATCH --error=./slurmoutputs/woltkaogu.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

# Load woltka
# Latest version on UGA GACRC Sapelo2 cluster is

module purge
module load woltka/0.1.5-QIIME2-2023.7
module list

# Define the input directory for sam files that will be used
mkdir ./alignments/${indir}/samfiles
mkdir ./woltka

# Woltka will use reads that were mapped concordantly and discordantly but I do not want that- I want to be very sure that all reads are bacterial and match...

mv ./alignments/$indir/*_bothreadsmapped.sam ./alignments/$indir/samfiles/

# The way that this has been installed on the GACRC opening this module allows you to call the standalone woltka program 'woltka', and also to use woltka through its qiime2 plug-in with 'qiime woltka'

woltka classify -i ./alignments/${indir}/samfiles \
-o ./woltka/ogu.biom

#-----------------------------WoltkaOGU.sh END----------------------------------

# We have now derived OGUs from our sequences which we can now assign taxonomy to using the greengenes2 database
# This feels like an appropriate choice because the gg2 taxonomy draws from the expansive GTDB, and can be used to aid in aligning shotgun data to amplicon data

# This is one of our few scripts where there's no need to use a bunch of fun dummy variables, theres really only one place this data source can be coming from

# Submit the taxonomy assignment script
sbatch ./scripts/ogutaxonomygg2.sh

#---------------------------ogutaxonomygg2.sh START-----------------------------

#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=ogutaxonomygg2
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=100G
#SBATCH --output=./slurmOutputs/ogutaxonomygg2.out
#SBATCH --error=./slurmOutputs/ogutaxonomygg2.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24

# load required modules
# module versions available on the UGA GACRC Server updated OCT 2023
module purge
module load QIIME2/2023.7
module load q2-greengenes2/2023.3-QIIME2-2023.7
module list

# download the greengenes2 taxonomy data
mkdir -p ./databases/greengenes2
wget http://ftp.microbio.me/greengenes_release/current/2024.09.taxonomy.asv.nwk.qza
mv 2024.09.taxonomy.asv.nwk.qza databases/greengenes2

# import woltka biom output file
qiime tools import \
	--input-path ./woltka/ogu.biom \
	--output-path ./woltka/ogu.biom.qza \
	--type FeatureTable[Frequency]

# filter out WoLr2 features not present in Greengenes2 and pull taxonomy from the phylogeny labels
#filter ogu table
qiime greengenes2 filter-features \
     --i-feature-table ./woltka/ogu.biom.qza \
     --i-reference ./databases/greengenes2/2024.09.taxonomy.asv.nwk.qza \
     --o-filtered-feature-table ./woltka/ogu.biom.gg2.filt.qza
# get taxonomy from labels
qiime greengenes2 taxonomy-from-table \
     --i-reference-taxonomy ./databases/greengenes2/2024.09.taxonomy.asv.nwk.qza \
     --i-table ./woltka/ogu.biom.gg2.filt.qza \
     --o-classification ./woltka/ogu.biom.gg2.filt.taxonomy.qza

# export taxonomy file
qiime tools export --input-path ./woltka/ogu.biom.gg2.filt.taxonomy.qza --output-path ./woltka/
# Rename the exported taxonomy file
mv ./woltka/taxonomy.tsv ./woltka/ogutaxonomy.gg2.tsv
# Re-format the taxonomy file in an actually sensible way...
# Define a header line for the re-formatted file
echo -e "Feature ID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tConfidence" > ./woltka/finalogutax.gg2.tsv
# Process each line of the file, skipping the header
tail -n +2 ./woltka/ogutaxonomy.gg2.tsv | while IFS=$'\t' read -r feature_id taxon confidence; do
    # Split the taxon field by "; " and remove the prefixes (d__, p__, etc.)
    domain=$(echo "$taxon" | awk -F'; ' '{print $1}' | sed 's/^d__//')
    phylum=$(echo "$taxon" | awk -F'; ' '{print $2}' | sed 's/^p__//')
    class=$(echo "$taxon" | awk -F'; ' '{print $3}' | sed 's/^c__//')
    order=$(echo "$taxon" | awk -F'; ' '{print $4}' | sed 's/^o__//')
    family=$(echo "$taxon" | awk -F'; ' '{print $5}' | sed 's/^f__//')
    genus=$(echo "$taxon" | awk -F'; ' '{print $6}' | sed 's/^g__//')
    species=$(echo "$taxon" | awk -F'; ' '{print $7}' | sed 's/^s__//')

    # Append the processed line to the final output file
    echo -e "$feature_id\t$domain\t$phylum\t$class\t$order\t$family\t$genus\t$species\t$confidence" >> ./woltka/finalogutax.gg2.tsv
done

# Convert the .tsv to .csv
tr '\t' ',' < ./woltka/finalogutax.gg2.tsv > ./woltka/finalogutax.gg2.csv

# export ogu table file
qiime tools export --input-path ./woltka/ogu.biom.gg2.filt.qza --output-path ./woltka/
mv ./woltka/feature-table.biom ./woltka/finalogutab.gg2.biom

# Clean up modules... biom-format and qiime don't play nicely on our hpc
module purge
ml biom-format/2.1.14-foss-2022a
module list

# convert ogu table file from .biom to .tsv
biom convert -i ./woltka/finalogutab.gg2.biom -o ./woltka/finalogutab.gg2.tsv --to-tsv
sed -i '1d' ./woltka/finalogutab.gg2.tsv
tr '\t' ',' < ./woltka/finalogutab.gg2.tsv > ./woltka/finalogutab.gg2.csv

# Now we'll count sample size across our final dataset
# Create a temporary file for the modified TSV without the first column
temp_tsv=$(mktemp)
# Remove the first column and save to the temporary file
cut -f2- ./woltka/finalogutab.gg2.tsv > "$temp_tsv"
# Extract the new header (sample IDs) from the modified TSV
header=$(head -n 1 "$temp_tsv")
IFS=$'\t' read -r -a sample_ids <<< "$header"
# Initialize the output CSV with a header
echo "SampleID,Sum" > countsummary_woltkaogu.csv
# Loop over each column index in the modified TSV
for i in $(seq 1 ${#sample_ids[@]}); do
  # Sum the numerical values in the column (skip the header)
  sum=$(awk -F'\t' -v col="$i" 'NR > 1 {sum += $col} END {print sum}' "$temp_tsv")
  # Add the sample ID and sum to the output CSV
  echo "${sample_ids[$i-1]},$sum" >> countsummary_woltkaogu.csv
done
# Clean up temporary file
rm "$temp_tsv"

#-------------------------greengenes2taxonomy.sh END----------------------------


# ************************** Kraken2 Taxonomy **********************************

# Download the required databases for kraken2
# Bacteria
sbatch --export=library=bacteria,dbname=microbesinclusive scripts/kraken2dl.sh
# Archaea
sbatch --export=library=archaea,dbname=microbesinclusive scripts/kraken2dl.sh
# Fungi
sbatch --export=library=fungi,dbname=microbesinclusive scripts/kraken2dl.sh
# UniVec_Core
sbatch --export=library=UniVec_Core,dbname=microbesinclusive scripts/kraken2dl.sh

#----------------------------kraken2dl.sh START---------------------------------
#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=kraken2dl
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=100G
#SBATCH --output=./slurmoutputs/kraken2dl.out
#SBATCH --error=./slurmoutputs/kraken2dl.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24

# Load kraken2 program
module purge
module load Kraken2/2.1.3-gompi-2022a

mkdir -p databases
mkdir -p databases/kraken2

kraken2-build --download-library $library --db ./databases/kraken2/$dbname

#----------------------------kraken2dl.sh END-----------------------------------

#----------------------------kraken2build.sh START------------------------------
#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=kraken2dl
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=500G
#SBATCH --output=./slurmoutputs/kraken2build.out
#SBATCH --error=./slurmoutputs/kraken2build.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

# Load kraken2 program
module purge
module load Kraken2/2.1.3-gompi-2022a

kraken2-build --build --db ./databases/kraken2/$dbname --clean --threads 128

#----------------------------kraken2build.sh END--------------------------------

sbatch ./scripts/mergecountsummary.sh


#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=mergecountsummary
#SBATCH --ntasks=1
#SBATCH --time=30:00
#SBATCH --mem=20G
#SBATCH --output=./slurmOutputs/mergecountsummary.out
#SBATCH --error=./slurmOutputs/mergecountsummary.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12

#!/bin/bash

# Output file
output="merged_output.csv"

# Temporary directory for processing
temp_dir=$(mktemp -d)

# Process each file and handle headers and sample counts
for file in countsummary_*.csv; do
    # Extract unique identifier from filename
    identifier=$(basename "$file" .csv | cut -d'_' -f2)

    # Extract the header (Sample Name, Read Count)
    head -n 1 "$file" > "$temp_dir/header_$identifier.csv"

    # Extract the 'Sample Name' and 'Read Count' columns, skipping header
    tail -n +2 "$file" | cut -d, -f1,2 > "$temp_dir/data_$identifier.csv"

    # Rename the second column (Read Count) with the unique identifier
    sed -i "1s/.*/Sample Name,$identifier/" "$temp_dir/data_$identifier.csv"
done

# Combine all headers (we need to do this separately)
header_files=($temp_dir/header_*.csv)
paste -d, "${header_files[@]}" | head -n 1 > "$output"

# Combine all data files by aligning columns
data_files=($temp_dir/data_*.csv)
paste -d, "${data_files[@]}" >> "$output"

# Clean up temporary directory
rm -r "$temp_dir"













###################################################3########
############################################################
################ 16S rDNA Amplicon Pipeline ################
############################################################
############################################################
############################################################

#******************************************************#
#***************** Initial Processing *****************#
#******************************************************#

# trim primer sequences from samples in parallel

# IMPORTANT: Update 'nSamples' to match your number of samples for parallelization to be split across.

nSamples=30

# Submit cutadapt script to trim primers

sbatch --array=1-$nSamples%30 seqtrim_515F806R.sh

-------SLURM SCRIPT seqtrim_515F806R.sh START-------------
#!/bin/bash
#SBATCH
#SBATCH --array=1-2%30
#SBATCH --partition=iob_p
#SBATCH --job-name=seqTrim
#SBATCH --ntasks=8
#SBATCH --time=4:00:00
#SBATCH --mem=40GB
#SBATCH --output=./slurmOutputs/seqTrim16S.%a.out
#SBATCH --error=./slurmOutputs/seqTrim16S.%a.out

# Script to trim primer sequences
# Assumes we are in the working directory
# Assumes working directory has the necessary subdirectories for storing output:
# ./trimmed_16S/ # to hold trimmed sequence files
# ./slurmOutputs/ # to hold output

# Get input files (paired forward/reverse reads)

FWD=`ls ./raw_16S/ | grep "_R1_" | head -n $SLURM_ARRAY_TASK_ID | tail -1` # get one forward read fi$
REV=`echo $FWD | sed 's/_R1_/_R2_/'` # get paired reverse read

# Extract sample name from the full input filename:
SAMPLE=$(echo ${FWD} | sed "s/_R1_\001\.fastq.gz*//")

#Load cutadapt
#Current version (OCT 2023) on GACRC is 4.5
module load cutadapt/4.5-GCCcore-11.3.0

# Trim primers
### Forward primer = 515F (Parada) = TGTGYCAGCMGCCGCGGTAA
### Reverse primer = 806R (Apprill) = GGACTACNVGGGTWTCTAAT
### Input comes from ./raw/ subdirectory
### Output goes into ./trimmed/ subdirectory
### Reads with no primer matches get discarded

#Cutadapt script
cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed --output=./trimmed_16S/${SAMPLE}_R1_trimmed.fastq.gz --paired-output=./trimmed_16S/${SAMPLE}_R2_trimmed.fastq.gz ./raw_16S/$FWD ./raw_16S/$REV

--------SLURM SCRIPT seqtrim_515F806R.sh END--------------

# Perform fastQC quality checking

sbatch fastQC.sh

------------SLURM SCRIPT fastQC.sh START------------------
#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=fastqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --mem=250G
#SBATCH --output=./slurmOutputs/fastQC.out
#SBATCH --error=./slurmOutputs/fastQC.out

# 1. Fire up fastQC:
module load FastQC/0.11.9-Java-11 #SEP 2023 current version on Sapelo2

# 2. Generate fastQC summaries for R1 and R2:
zcat ./trimmed_16S/*R1_trimmed.fastq.gz | fastqc -t 16 -o ./fastQCoutput_16S stdin:R1
zcat ./trimmed_16S/*R2_trimmed.fastq.gz | fastqc -t 16 -o ./fastQCoutput_16S stdin:R2

-------------SLURM SCRIPT fastQC.sh END-------------------


# download fastQC outputs (./fastQCoutput/*.html)
# './' points to your current working directory on your home machine

scp pab14613@xfer.gacrc.uga.edu:/scratch/pab14613/Y1P/fastQCoutput/R1_fastqc.html ./ #R1
scp pab14613@xfer.gacrc.uga.edu:/scratch/pab14613/Y1P/fastQCoutput/R2_fastqc.html ./ #R2

# view quality metrics on a web browser by opening the HTML file
# you can use the quality drop-off to determine where you will trim your sequences for later steps

#******************************************************#
#**************** ASV / OTU Generation ****************#
#******************************************************#

#************ DADA2 Processing and Cleanup ************#

# Carry out ASV denoising using the R implementation of dada2, and export a biom file for QIIME2 import.

sbatch --export=R1cutoff=250,R2cutoff=250 trimmedtoASVs.sh # cutoff = 0 means no trimming

#******** 97% OTU Generation and Sample Cleanup *******#

# Use QIIME2/vsearch to perform OTU aggregation, carry out additional UCHIME chimera detection, and filter low abundance OTUs and low depth samples.

sbatch ASVs2OTUs.sh

-------SLURM SCRIPT ASVs2OTUs.sh START-------------

-------SLURM SCRIPT greengenes2amplicontaxonomy.sh START-------------

#******************************************************#
#************** Taxonomic Classification **************#
#******************************************************#

# we will assign taxonomy to the OTUs based on  phylogenetic placement through greengenes2 so that our taxonomy database matches that which we use for the bacterial assignments.

# Following taxonomic assignment we will remove any non-target sequences (e.g., mitochondria, chloroplast, non-bacterial).

sbatch greengenes2amplicontaxonomy.sh

-------SLURM SCRIPT greengenes2amplicontaxonomy.sh START-------------
#!/bin/bash
#SBATCH
#SBATCH --partition=batch
#SBATCH --job-name=gg2amplicontaxonomy
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=250GB
#SBATCH --output=./slurmOutputs/gg2amplicontaxonomy.%a.out
#SBATCH --error=./slurmOutputs/gg2amplicontaxonomy.%a.out

# download full length 16S rRNA sequence backbone phylogeny

wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza

# map your sequences to the greengenes2 backbone phylogeny

qiime greengenes2 non-v4-16s \
	-i-table <BIOMFILE.biom.qza> \
	-i-sequences <REPSETFILE.fna.qza> \
	-i-backbone 2022.10.backbone.full-length.fna.qza \
	-o-mapped-table <OUTPUTBIOMFILE.gg2.biom.qza> \
	-o-representatives <OUTPUTREPSET.gg2.fna.qza>

# classify taxonomy of mapped sequences

qiime greengenes2 taxonomy-from-table \
	--i-reference-taxonomy 2022.10.taxonomy.asv.nwk.qza \
	--i-table <OUTPUTBIOMFILE.gg2.biom.qza> \
	--o-classification <OUTPUTTAXONOMY.gg2.taxonomy.qza>

# export required files for downstream analysis

-------SLURM SCRIPT greengenes2amplicontaxonomy.sh END-------------

###################################################3########
############################################################
################## Combined Data Analysis ##################
############################################################
############################################################
############################################################
