#!/usr/bin/env bash
#SBATCH --job-name="Adapter_removal"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mem=4G

#Data Preprocessing

#User needs to input path to this directory where the input data is (source) and the path to the directory where output data should go (destination) when running the script with sbatch
#example: sbatch adapter_removal.sh </path/to/source/> </path/to/destination1/> </path/to/destination2/>

# source/path where the input files are present.
source=$1 

#destination/path where the processed data should be stored.
destination1=$2

#destination/path where the report/result of the analysis sould be stored.
destination2=$3

#Load all the modules required for the QC analysis

module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Quality_control/cutadapt/2.5
module add UHTS/Analysis/MultiQC/1.8


mkdir ${destination1}/trimmed_data

cd ${destination1}/trimmed_data

#Create a link to the raw files for data processing 
#for loop that loops over all the files with .fsatq extension in the given path and creates individual links for all the .fastq files in the current directory

for file in ${source}/*.fastq.gz; do ln -s "$file" . ; done

# Code for Clipping 3' adapter
# 1> is used to send the standard output of the code to the filename that follows

for i in `ls -1 *.fastq.gz`;
do cutadapt -j 0 -a CTGTAGGCACCATCAAT -q 25 -m 25 --discard-untrimmed -o $(basename ${i} .fastq.gz)_clpd.fastq.gz \
${i} 1> $(basename ${i} .fastq.gz)_clpd_log.txt; rm $i; done

# for trimming 4 nt (randomized bases) from the 3' end

for i in `ls -1 *_clpd.fastq.gz`;
do cutadapt -j 0 -q 25 --cut -4 -m 25 -o $(basename ${i} .fastq.gz)_tr.fastq.gz ${i} 1> $(basename ${i} .fastq.gz)_tr_log.txt; done

#QC analysis of the processed reads

mkdir ${destination2}/QC_trimmed_reads
cd ${destination2}/QC_trimmed_reads

for j in ../*_tr.fastq.gz; do ln -s "$j" . ; done

for j in `ls -1 *_tr.fastq.gz`; do fastqc -t 32 $j; rm $j; done

multiqc .

#You can now download the html files on your local computer to assess the quality of the trimmed reads.

