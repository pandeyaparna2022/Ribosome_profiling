#!/usr/bin/env bash
#SBATCH --job-name="Mapping"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=10G

#User needs to input path to this directory where the input data is (source), the path to the directory where the indicies for allignment are and the path to the directory where output data should go (destination) when running the script with sbatch
#example: sbatch mapping.sh </path/to/source/> </path/to/indices/> </path/to/destination/>

source=$1
indices=$2
destination=$3

#Load all the modules required for the allignment
module load UHTS/Aligner/bowtie/1.2.0
module load UHTS/Analysis/samtools/1.10

# Mapping the trimmed reads to undesired RNAs and filtering them out

mkdir ${destination}/filtered_data
cd ${destination}/filtered_data/

for x in $(ls ${source}/*tr.fastq.gz); do gunzip ${x}; done
for x in ${source}*tr.fastq; \
do bowtie -t -p 4 ${indices}/GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb  ${x} --un $(basename ${x} .fastq)_no_r-t-sno-sn-RNA.fastq 2> $(basename ${x} .fastq.gz)_no_r-t-sno-sn-RNA_log.txt > /dev/null; done

for file in ${source}*tr.fastq; do gzip ${file}; done

#Mapping the filtered reads to the genome
mkdir ${destination}/mapped_bamfiles
cd ${destination}/mapped_bamfiles/

for x in ${destination}/filtered_data/*no_r-t-sno-sn-RNA.fastq; \
do bowtie -S -t -p 4 -v 1 -m 1 --best --strata ${indices}/GRCh38.p13.genome -q ${x} 2> $(basename ${x} .fastq)_GRCh38_p13_log.txt | samtools view -h -F 4 -b > $(basename ${x} .fastq)_GRCh38_p13.bam; done

# Sort the BAM file
for x in $(ls -d *.bam);do samtools sort -@ 4 ${x} -o $(basename ${x} .bam)_sorted.bam; done

# Remove the unsorted BAM file
rm *GRCh38_p13.bam


