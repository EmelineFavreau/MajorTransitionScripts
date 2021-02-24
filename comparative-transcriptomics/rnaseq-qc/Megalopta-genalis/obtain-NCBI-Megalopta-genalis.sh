#!/bin/bash

# Obtain Megalopta genalis raw input from NCBI
# Copyright Emeline Favreau, UCL



# download the data BioProject-331103-SraAccList.txt
~/programs/sratoolkit.2.10.8-centos_linux64/bin/prefetch --option-file BioProject-331103-SraAccList.txt

ls -lhrt /lustre/scratch/scratch/ucfaeef/ncbi/sra/

# convert to fastq 
~/programs/sratoolkit.2.10.8-centos_linux64/bin/fastq-dump \
	--outdir /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/${Genusspecies}/input/ \
	--split-files /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948578.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948576.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948574.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948571.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948569.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948567.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948565.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948561.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948522.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948523.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948524.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948525.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948530.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948532.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948534.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948535.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948538.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948541.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948543.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948545.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948546.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948549.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948551.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948553.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948555.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948557.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948559.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948573.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948580.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR3948582.sra


ls -lrth /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/${Genusspecies}/input/

# compress them all 
gzip /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/${Genusspecies}/input/*.fastq

cd inputOnScratch

# get the genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/865/705/GCF_011865705.1_USU_MGEN_1.2/GCF_011865705.1_USU_MGEN_1.2_genomic.fna.gz

# get the gff
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/865/705/GCF_011865705.1_USU_MGEN_1.2/GCF_011865705.1_USU_MGEN_1.2_genomic.gff.gz