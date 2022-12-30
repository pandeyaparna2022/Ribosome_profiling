#!/usr/bin/env bash
#SBATCH --job-name="QC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mem=8G


#QC analysis

#User needs to input path to this directory where the input data is (source) and the path to the directory where output data shuld go (destination) when running the script with sbatch
#example: sbatch QC.sh </path/to/source/> </path/to/destination>

source=$1
destination=$2

#load required modules
module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Analysis/MultiQC/1.8

#change directory to interim
cd ${destination}


#Create a link to the raw files, in the current directory for qc analysis
#code: for loop that loops over all the files with .fsatq.gz extension in the given path and creates individual links for all the .fastq.gz files in the current directory

for file in ${source}*.fastq.gz; do ln -s "$file" . ; done

for i in `ls -1 *.fastq.gz`;
do fastqc -t 6 $i; rm $i; 
done

multiqc .


#You can now download the html files on your local computer to assess the quality of raw reads
