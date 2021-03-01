#!/bin/bash

# Obtain longest transcript for  each gene in Ceratina calcarata
# Copyright Emeline Favreau, UCL

## Analysis overview
### Step 1: get the gene and protein ids from gff files
### Step 2: get the protein length from the proteomes
### Step 3: obtain the longest isoform for each protein
### Step 4: subset protein fasta file for just those isoforms

#################################################################################
## Step 1: get the gene and protein ids from gff files

# obtain gff
ssh ucfaeef@transfer02

cd /lustre/scratch/scratch/ucfaeef/orthology-analysis/input/gff

cp /mnt/gpfs/live/rd01__/ritd-ag-project-rd011y-eefav92/Ceratina_calcarata/ShellRehan2019/gff/GCF_001652005.1_ASM165200v1_genomic.gff.gz .

mv GCF_001652005.1_ASM165200v1_genomic.gff.gz Ceratina_calcarata.gff.gz


exit 

thisspecies="Ceratina_calcarata"

cd ../../

# check that gff is NCBI format (RefSeq)
# zcat inputOnScratch/gff/${thisspecies}.gff.gz | head
# here maybe I need to include only gene from third column
zcat inputOnScratch/gff/${thisspecies}.gff.gz | grep "CDS" | head
head inputOnScratch/${thisspecies}.faa

# this is specific to each GFF (not a norm unless from Ensembl)
# get gene and protein IDs from GFF file, e.g.108630137	XP_017888733.1
# use uniq to remove duplicated
# useful link : https://www.biostars.org/p/388936/
# remove in column 9 anything that is before and includes ;Dbxref=GeneID:
# remove in column 9 anything that is after ;Name
# remove the string ,Genbank: and replace it by a tab
zcat inputOnScratch/gff/${thisspecies}.gff.gz | awk '$3=="CDS"{gsub(".+;Dbxref=GeneID:", "", $9); gsub(";Name=.+", "", $9); gsub(",Genbank:", "\t", $9); print $9}' | sort | uniq > inputOnScratch/gff/${thisspecies}.tsv


#################################################################################
## Step 2: get the protein length from the proteomes

# first, tidy fastaa headers by checking first ( head -n 1 inputOnScratch/${thisspecies}.faa | cut -d " " -f 1)
# remove the second part of the sequence header after the gene ID, the delimmiter is space
# fx2tab calculate the sequence length
# sort the result
~/programs/./seqkit fx2tab -n -l inputOnScratch/${thisspecies}.faa \
   | awk 'BEGIN{FS="\t"}{gsub(" .*", "", $1); print}' \
   | sort \
  > inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv



#################################################################################
## Step 3: obtain the longest isoform for each protein
# based on Anindita's longest_protein_isoform_Step_2.Rmd
module unload compilers
module unload mpi
module load r/recommended

Rscript --vanilla subset-longest-protein-isoform.R inputOnScratch/gff/${thisspecies}.tsv inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt



#################################################################################
## Step 4: subset protein fasta file for just those isoforms
# subset protein sequences for those longest proteins

# update the faa headers for short name (ie just the gene name)
awk '{gsub(" .*", ""); print}' inputOnScratch/${thisspecies}.faa > inputOnScratch/${thisspecies}-short-name.faa

# remove protein sequences that are not the longest isoform
seqtk subseq inputOnScratch/${thisspecies}-short-name.faa \
	inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt \
	> inputOnScratch/primary_transcripts/${thisspecies}-longest-isoforms.faa

# this protein fasta has the longest isoform, with headers containing just gene id (i.e. >XP_014597714.1)