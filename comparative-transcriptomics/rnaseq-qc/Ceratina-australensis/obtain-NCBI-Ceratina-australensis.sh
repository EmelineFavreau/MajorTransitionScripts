#!/bin/bash

# Obtain Ceratina australensis raw input from NCBI
# Copyright Emeline Favreau, UCL



# download the data bioproject-PRJNA302037-accession-list
~/programs/sratoolkit.2.10.8-centos_linux64/bin/prefetch --option-file bioproject-PRJNA302037-accession-list

ls -lhrt /lustre/scratch/scratch/ucfaeef/ncbi/sra/

~/programs/sratoolkit.2.10.8-centos_linux64/bin/fastq-dump \
	--outdir /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/Ceratina-australensis/input/ \
	--split-files /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2917198.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2917152.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916699.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916653.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916652.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916648.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916647.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916637.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916631.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916026.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2916025.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR2915407.sra 


ls -lrth /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/Ceratina-australensis/input/

# compress them all 
gzip /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/Ceratina-australensis/input/*.fastq

# we use the genome from Sandra Rehan
Ceratina australensis genome: 
https://drive.google.com/drive/folders/0B1P0qnA-7O3Afm1ieUZKa2VhSGpacWVsNXVWWGZ0eVNJLUFfX0xZTG9YRWs1YUhoZS1qb00?usp=sharing

gzip australensis_DNA.noNUL.fasta

# we use the annotation from Sandra Rehan
https://drive.google.com/drive/folders/0B1P0qnA-7O3Afm1ieUZKa2VhSGpacWVsNXVWWGZ0eVNJLUFfX0xZTG9YRWs1YUhoZS1qb00?usp=sharing

gzip australensis.MAKER_run_v1.cfig_v5.MAKER_ONLY_FORMATTED.gff.filtered