#Course work
Course: RNA- sequencing course, module: Ribosome_profiling


These codes are written for analysis of Ribo-Seq data.

Recommendation:

Prior to running the code it is suggested to have the following structured layout for the data and results. This is not essential but if you deviate from the structure you might need to make changes to the paths within the scripts wherever necessary

Create 2 folders data and reports.

All data including Raw fasta files, intermediate processed files and annotation files or precursors to annotation files are kept withing ordered subfolders within data folder.

../data/raw - raw fastq.gz files

../data/external -all the data files that are downloaded from various outside sources used for annotation.

   /annotation_genome - files for annotating the genome and building genomic indices
  
    /indices - genomic indices
    
  /annnotation_transcriptome - files for annotating the transcriptome and building transcriptomic indices 
  
    /indices - transcriptomic indices
    
  /annotation_undesired_RNA - files for annotating the undesired RNA and building indices for undesired RNAs
  
    /indices - indices for undesired RNAs
    
  
All results and reports , intermediate and final are directed to the ordered folders within the reports folder.





