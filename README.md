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
All the scripts require at least 1 input path and 1 output path indicating the location of the input and the output file respectively. The script can be run as follows:

example.sh <path/to/input/files/> <path/to/output/files/>


In caseses where different types of files are required as input two input paths are required. The script can be run as follows:

example.sh <path/to/type1_input/files/> <path/to/type2_input/files/> <path/to/output/files/>








