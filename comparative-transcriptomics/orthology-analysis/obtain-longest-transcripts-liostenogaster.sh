#!/bin/bash

# Obtain longest transcript for  each gene in Liostenogaster
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

cd /home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/orthology-analysis/inputOnScratch/gff

# GFF has 170,079 features
cp /mnt/gpfs/live/ritd-ag-project-rd00qm-mbent70/Liostenogaster_flavolineata/genome/L.flavolineata.evm.consensus.annotation.v2a.gff3.uni .

exit

mv inputOnScratch/gff/L.flavolineata.evm.consensus.annotation.v2a.gff3.uni inputOnScratch/gff/Liostenogaster_flavolineata.gff

gzip inputOnScratch/gff/Liostenogaster_flavolineata.gff

thisspecies="Liostenogaster_flavolineata"

# check that gff is NCBI format (RefSeq) here it is not
# zcat ${thisspecies}.gff.gz | head

# this is specific to each GFF (not a norm unless from Ensembl)
# get gene and protein IDs from GFF file, e.g. Lflavo2a009244P1 and Lflavo2a009244T1_cds 
# use uniq to remove duplicated

# Lflavo2a009244P1;Target=Lflavo2a009244P1 1 3;ID=Lflavo2a009244T1_cds;Name=Lflavo2a009244T1_cds
# useful link : https://www.biostars.org/p/388936/
# remove in column 9 Parent=
# remove in column 9 everything between target and name, but print three columns
# remove in this output everything between target and name (this time, the spaces became underscore)
# remove partial
zcat inputOnScratch/gff/${thisspecies}.gff.gz | awk '$3=="CDS"{gsub("Parent=", "", $9); gsub(";Target=.*;Name=", " ", $9) ; print $9 $10 $11}' | awk '{gsub(";Target=.*;Name=", "\t"); print}' | grep -v "partial" | sort | uniq > inputOnScratch/gff/${thisspecies}.tsv

# number of proteins (this includes potential multiple proteins for each gene)
# 13914


#################################################################################
## Step 2: get the protein length from the proteomes

# in inputOnScratch

# fx2tab calculate the sequence length for each protein
# e.g. Lflavo2a000001P1	protein is 429 amino acids long
# sort the result
~/programs/./seqkit fx2tab -n -l inputOnScratch/${thisspecies}.faa \
  | sort \
  > inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv

# 17208 protein sequences

# For unknown reason to me,the fasta file has more sequences than the gff
# so I need to intersect gff and fasta.

# get names of protein from gff file (13,914)
cut -f 1 inputOnScratch/gff/${thisspecies}.tsv > tmp/${thisspecies}-gff-names

# get names of protein from fasta file (17,208)
cut -f 1 inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv | sort > tmp/${thisspecies}-fasta-names

# get common protein names (13,914)
comm -12 tmp/${thisspecies}-gff-names tmp/${thisspecies}-fasta-names > tmp/${thisspecies}-common-names

# get protein unique to fasta: (3,294)
comm -13 tmp/${thisspecies}-gff-names tmp/${thisspecies}-fasta-names > tmp/${thisspecies}-names-unique-to-fasta


# reduce the file with protein names and length to keep only the common proteins (13,898)
grep -v -f tmp/${thisspecies}-names-unique-to-fasta inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv > inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv

# keep only the proteins that are common to both
grep -f tmp/${thisspecies}-common-names inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv > inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered-matching-gff.tsv

# check number of unique protein names: 13898
cut -f1 inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered-matching-gff.tsv > tmp/protein-names-13898

# keep proteins in gff that are also found in fasta
grep -f tmp/protein-names-13898 inputOnScratch/gff/${thisspecies}.tsv > inputOnScratch/gff/${thisspecies}-filtered.tsv

# this filtering leaves fewer lines: 13898
# wc -l inputOnScratch/gff/${thisspecies}-filtered.tsv
# wc -l inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv


# how much length did this remove? from 8,977,259 to 6,287,089 (a third)
# awk '{SUM+=$2}END{print SUM}' inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv
# awk '{SUM+=$2}END{print SUM}' inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv


#################################################################################
## Step 3: obtain the longest isoform for each protein
# based on Anindita's longest_protein_isoform_Step_2.Rmd
module unload compilers
module unload mpi
module load r/recommended

Rscript --vanilla subset-longest-protein-isoform.R inputOnScratch/gff/${thisspecies}-filtered.tsv inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt


#################################################################################
## Step 4: subset protein fasta file for just those isoforms
# subset protein sequences for those longest proteins

# remove protein sequences that are not the longest isoform
seqtk subseq inputOnScratch/${thisspecies}.faa \
	inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt \
	> inputOnScratch/primary_transcripts/${thisspecies}-longest-isoforms.faa

# this protein fasta has the longest isoform, with headers containing just gene id (i.e. >Lflavo2a000001P1)