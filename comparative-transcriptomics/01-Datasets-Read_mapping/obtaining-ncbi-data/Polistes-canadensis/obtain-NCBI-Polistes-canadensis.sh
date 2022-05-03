# Obtain Polistes canadensis Patalano raw input from NCBI
# Copyright Emeline Favreau, UCL

Genusspecies="Polistes-canadensis"
dataset="Patalano"

# download the data SraAccList.txt
~/programs/sratoolkit.2.10.8-centos_linux64/bin/prefetch --option-file SraAccList.txt

ls -lhrt /lustre/scratch/scratch/ucfaeef/ncbi/sra/

# convert to fastq 
~/programs/sratoolkit.2.10.8-centos_linux64/bin/fastq-dump \
	--outdir /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/${Genusspecies}/input/${dataset}/ \
	--split-files /lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519108.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519109.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519110.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519111.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519112.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519113.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519114.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519115.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519116.sra \
	/lustre/scratch/scratch/ucfaeef/ncbi/sra/SRR1519117.sra

#ls -lrth /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/${Genusspecies}/input/${dataset}

# compress them all 
gzip /lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/${Genusspecies}/input/${dataset}/*.fastq

cd inputOnScratch

# get the genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/313/835/GCF_001313835.1_ASM131383v1/GCF_001313835.1_ASM131383v1_genomic.fna.gz

# get the gff
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/313/835/GCF_001313835.1_ASM131383v1/GCF_001313835.1_ASM131383v1_genomic.gff.gz

