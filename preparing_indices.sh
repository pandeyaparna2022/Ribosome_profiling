#!/usr/bin/env bash
#SBATCH --job-name="Preparing_Indices"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --mem=10G

# load bowtie 
module load UHTS/Aligner/bowtie/1.2.0
module load SequenceAnalysis/blat/36

#User needs to input path to this directory where the genomic files are (source) (data/external). The prepared indicies will be stored in a new folder within the source folder called indices
#example: sbatch preparing_Indices.sh </path/to/source/> 

source=$1
destination=$2

cd ${source}

#Use bowtie build to make indices

# For the genome
cd annotation_genome/
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bowtie-build Homo_sapiens.GRCh38.dna.primary_assembly.fa ${destination}/GRCh38.p13.genome
faToTwoBit Homo_sapiens.GRCh38.dna.primary_assembly.fa  GRCh38.p13.genome.2bit # done once for Riboseq QC analysis (after mapping to the genome)
gzip Homo_sapiens.GRCh38.dna.primary_assembly.fa

cd ${source}

# For the "undesired" RNAs
cd annotation_undesired_RNA/

#Combine all the different files containing the different types of undesired RNA. This is done only once.
cat *.txt > GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb.fa 

bowtie-build GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb.fa ${destination}/GRCh38_p13_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb

cd ${source}

# For the transcriptome
cd annotation_transcriptome/
bowtie-build GRCh38_p13_APPRIS_CDS_plus18.fa ${destination}/GRCh38_p13_APPRIS_CDS_plus18


#Converthing into single-line format for generating codon occupancy plots. This is done only once.

awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < GRCh38_p13_APPRIS_CDS_plus18.fa > GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa

cp GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa ${destination}/
