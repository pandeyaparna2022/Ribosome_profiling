#Course work
Course: RNA- sequencing course, module: Ribosome_profiling

These codes are written for analysis of Ribo-Seq data. Please make sure you have access to the raw sequencing data and have downloaded the individual annotation files prior to data analysis. A step-by-step guide for downloading the annotation files are provided at file:///C:/Users/pandapar/Desktop/Uni_courses_2022/RNA_sequencing/RNASeq_2022_RP.html#prepare-annotations under Prepare annotations section (author: Dr. Puneet Sharma)

Recommendation:

Prior to running the code it is suggested to have the following structured layout for the data and results. This is not essential, just a recommentation.

Create a folder named data.

All data including Raw fasta files, intermediate processed files and annotation files or precursors to annotation files are kept in ordered subfolders within data folder.

../data/raw - for raw fastq.gz files

../data/external -for all the data files used for annotation that are downloaded from various web resources.
../data/external/annotation_genome - for files for annotating the genome and building genomic indices
../data/external/annnotation_transcriptome -for  files for annotating the transcriptome and building transcriptomic indices 
../data/external/annotation_undesired_RNA -for  files for annotating the undesired RNA and building indices for undesired RNAs
../data/external/indices - for indices that are made from the downloaded annnotations
   
    
../data/interim - for all the processed data
 
 
All the scripts are written with SBATCH options at the beginning so that the script could be submitted to the SLURM workload manager to reduce the conputational load at the head node of the IUB cluster.

All the scripts require at least 1 input path and 1 output path indicating the location of the input file to be processed and the location of the output file respectively. The script can be run as follows:

example.sh <path/to/input/files/> <path/to/output/files/>

Following scripts require 1 input and 1 output path:

QC.sh </path/to/input/fastq.gz/files/> </path/to/destination/for/output/files>
adapter_trimming.sh </path/to/input/fastq.gz/files/> </path/to/destination/for/output/files>
preparing_indices.sh </path/to/annotation/files/> </path/to/destination/for/output/files/>


In some caseses where two or more different types of input files are required, two input paths need to be given. For the scripts that require 2 input files the scripts can be run as follows:

example.sh <path/to/type1_input/files/> <path/to/type2_input/files/> <path/to/output/files/>

Following scripts require 2 input and 1 output paths:
 
mapping.sh </path/to/input/trimmed/fastq.gz/files/> <path/to/prepared/indices> </path/to/output/files/>
feature_counts.sh </path/to/annotation/file(.gtf)/for/transcriptome/> </path/to/sorted/bam/files/> </path/to/output/files/>

For the scripts that require more than 2 input files the scripts can be run as follows:
 
example.sh <path/to/type1_input/files/> <path/to/type2_input/files/> <path/to/type3_input/files/> <path/to/output/files/>
 
Following scripts require 3 input and 1 output paths:
Codon_occupancy_generation.sh <path/to/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa/file> <path/to/sam/files/generated/by mapping/filtered/reads/to/the/transcriptome/> <path/to/Codon_occupancy_cal.sh/script/> <path/to/output/files/>
 







