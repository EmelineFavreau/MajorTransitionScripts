#!/bin/bash

# Obtain Ceratina calcarata raw input from NCBI
# Copyright Emeline Favreau, UCL

cd inputOnScratch

# download the data 
~/programs/sratoolkit.2.10.8-centos_linux64/bin/prefetch --option-file bioproject-PRJNA434715-accession-list

# ls -lhrt /lustre/scratch/scratch/ucfaeef/ncbi/sra/

# convert to fastq 
~/programs/sratoolkit.2.10.8-centos_linux64/bin/fastq-dump --outdir /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/Ceratina-calcarata/input/ --split-files /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281997.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281999.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281996.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281983.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282001.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282002.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282006.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282000.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281992.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281990.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282008.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281981.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281993.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282011.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281991.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281986.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282005.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282009.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281989.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281985.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281987.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282010.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282004.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282003.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281984.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281979.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281988.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281995.sra  /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281982.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13282007.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281994.sra  /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281980.sra /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR13281998.sra


ls -lrth /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/Ceratina-calcarata/input/

# compress them all 
gzip /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/Ceratina-calcarata/input/*.fastq

# we use the genome from Sandra Rehan
# Ceratina calcarata genome: https://drive.google.com/drive/folders/0B1P0qnA-7O3AfmotYXQ1c2tKdzNTQ3FTc3plNkRvRkViaWcwRFFvV0dIREQ3bEVoM1JFQ1k?usp=sharing
#gzip Rehan1.Ccalc_v1.1.draft_build.over200.fa
# there is an official version now, so I will use this 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/652/005/GCF_001652005.1_ASM165200v1/GCF_001652005.1_ASM165200v1_genomic.fna.gz


# we use the annotation from Sandra Rehan
# Ceratina calcarata annotation: https://drive.google.com/drive/folders/0B1P0qnA-7O3AfmotYXQ1c2tKdzNTQ3FTc3plNkRvRkViaWcwRFFvV0dIREQ3bEVoM1JFQ1k?usp=sharing

#gzip CCalc_v1.1.draft_genome_remap.v5c_3G2.MAKERonly.gff.filtered
# there is an official version now, so I will use this 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/652/005/GCF_001652005.1_ASM165200v1/GCF_001652005.1_ASM165200v1_genomic.gff.gz

