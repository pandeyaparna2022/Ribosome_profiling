#!/usr/bin/env bash

#BATCH --job-name="Feature_counts"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=10G

module load UHTS/Analysis/subread/2.0.1

source1=$1 # directory where transcriptome annotation file is located
source2=$2 # directory where sorthed bam files are located
destination=$3 # directory where output should go

#change directory to source where the transcriptome  annotation file is located
cd ${source1}

#unzip the annotation file for transcriptome
gunzip Homo_sapiens.GRCh38.108.gtf.gz

#change directory to destination
cd ${destination}

#make a directory to store counts table and change to that directory
mkdir counts_table
cd counts_table

#use featurecounts to map the bam files to the coding sequence in the transcriptome annotation and collapse to yield one value per gene 
featureCounts -T 4 -t CDS -g gene_id -a ${source1}/Homo_sapiens.GRCh38.108.gtf -o CDS_counts_rawfile.txt ${source2}/*_GRCh38_p13_sorted.bam


#use featurecounts to map the bam files to the exonic sequence in the transcriptome annotation and collapse to yield one value per biotype 
featureCounts -T 4 -t exon -g gene_biotype -a ${source1}/Homo_sapiens.GRCh38.108.gtf -o biotype_counts_rawfile.txt ${source2}/*_GRCh38_p13_sorted.bam


#extract columns 1,7,8,9,10 from the raw counts table for further analyis
cut -f 1,7-10 CDS_counts_rawfile.txt > CDS_counts_processed.txt
cut -f 1,7-10 biotype_counts_rawfile.txt > biotype_counts_processed.txt

#change to source destination and compress the annotation file
cd ${source1}
gzip Homo_sapiens.GRCh38.108.gtf

