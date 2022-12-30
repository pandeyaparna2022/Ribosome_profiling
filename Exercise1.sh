#!/usr/bin/env bash
#SBATCH --job-name="QC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mem=8G

#Exercise_1
wd=$1 #working directory where you want to download the files and do further processing, user needs to input this directory when running the script with sbatch
# example sbatch Exercise1.sh </path/to/your/desired/working/directory/>

module load UHTS/Analysis/sratoolkit/2.10.7
module load UHTS/Quality_control/fastqc/0.11.9

cd ${wd}

# -- split-files is used to write reads into separate files for forward and reverse reads in case of paired end reads, this is not necessary for our data since this is single end reads.
# -- progress shows percentage of the download completed. Not essential.

fasterq-dump --split-files --progress SRR9596295
fasterq-dump --split-files --progress SRR9596296
fasterq-dump --split-files --progress SRR9596300
fasterq-dump --split-files --progress SRR9596303
fasterq-dump --split-files --progress SRR9596304
fasterq-dump --split-files --progress SRR9596310

#create a compressed file for any file ending in .fastq
gzip *.fastq

#Rename the files
mv ./SRR9596295.fastq.gz ./Biever_Somata_Poly_1.fastq.gz
mv ./SRR9596296.fastq.gz ./Biever_Somata_Poly_2.fastq.gz
mv ./SRR9596300.fastq.gz ./Biever_Somata_Poly_3.fastq.gz
mv ./SRR9596310.fastq.gz ./Biever_Neuropil_Poly_1.fastq.gz
mv ./SRR9596303.fastq.gz ./Biever_Neuropil_Poly_2.fastq.gz
mv ./SRR9596304.fastq.gz ./Biever_Neuropil_Poly_3.fastq.gz

#create and go to a directory (within the current directory) where you want to keep your QC results
mkdir QC_results_exercise1
cd QC_results_exercise1

#Create a link to the raw files (which are in your working directory), in the current directory (the one you just created) for qc analysis
#for loop that loops over all the files with .fsatq.gz extension in the given path and creates individual links for all the .fastq.gz files in the current directory

for file in ${wd}/*.fastq.gz; do ln -s "$file" . ; done

#for loop that loops over all the links to the files with .fsatq.gz extension in current directory and performs QC on the.fastq.gz files
for i in `ls -1 *.fastq.gz`;
do fastqc -t 6 $i; rm $i;
done

# Go back a directory where all the downloaded fastq.gz files are present
cd ..

# Deleate them as they take up a lot of space. We are deleting them here because we are not working with them anymore.
rm -r fastq.gz

# You can now download all the html files (using the wildcard *.html) onto your local computer and assess the quality of the reads
# example code : scp -r <username>@binfservms01.unibe.ch:/path/to/the/QC_results_exercise1/*.html /path/to/the/destination/folder/in/local/computer/
