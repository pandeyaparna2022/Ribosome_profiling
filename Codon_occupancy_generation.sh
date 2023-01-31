#!/usr/bin/env bash

#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=02:00:00
#SBATCH --job-name=Codon_occupancy_generation
#SBATCH --mail-user=aparna.pandey@students.unibe.cH
#SBATCH --mail-type=end,fail
#SBATCH --output=Codon_occupancy_generation_output_%j.txt
#SBATCH --error=Codon_occupancy_generation_error_%j.txt

source1=$1
source2=$2
source3=$3
destination=$4

cd ${destination}

${source3}/Codon_occupancy_cal.sh \
${source1}/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa \
${source2}/RPF_KO_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep1_Codon_occupancy.txt

${source3}/Codon_occupancy_cal.sh \
${source1}/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa \
${source2}/RPF_KO_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep2_Codon_occupancy.txt

${source3}/Codon_occupancy_cal.sh \
${source1}/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa  \
${source2}/RPF_WT_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep1_Codon_occupancy.txt

${source3}/Codon_occupancy_cal.sh \
${source1}/GRCh38_p13_APPRIS_CDS_plus18_SingleLine.fa \
${source2}/RPF_WT_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p13_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep2_Codon_occupancy.txt
