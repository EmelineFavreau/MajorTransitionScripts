#!/bin/bash

# Obtain longest transcript for  each gene in Ceratina australensis
# Copyright Emeline Favreau, UCL

## Analysis overview
### Step 1: get the gene and protein ids from gff files
### Step 2: get the protein length from the proteomes
### Step 3: obtain the longest isoform for each protein
### Step 4: subset protein fasta file for just those isoforms

#################################################################################
## Step 1: get the gene and protein ids from gff files


thisspecies="Ceratina_australensis"

# this is specific to each GFF (not a norm unless from Ensembl)
# get gene and protein IDs from GFF file, e.g. Caust.v2_021703-RA       Caust.v2_021703-RA
# use uniq to remove duplicated
# useful link : https://www.biostars.org/p/388936/
# remove in column 9 anything that is before and includes ;Dbxref=GeneID:
# remove in column 9 anything that is after ;Name
# remove the string ,Genbank: and replace it by a tab
zcat inputOnScratch/gff/${thisspecies}.gff.gz | awk '$3=="gene"{gsub("ID=", "", $9); gsub(";$", "", $9); gsub(";Name=.*;Alia
s=", "\t", $9); print $9}' | sort | uniq > inputOnScratch/gff/${thisspecies}.tsv


#################################################################################
## Step 2: get the protein length from the proteomes

# first, tidy fastaa headers by checking first ( head -n 1 inputOnScratch/${thisspecies}.faa | cut -d " " -f 1)
# remove the second part of the sequence header after the gene ID, the delimmiter is space
# fx2tab calculate the sequence length
# sort the result
~/programs/./seqkit fx2tab -n -l inputOnScratch/${thisspecies}.faa \
   | awk 'BEGIN{FS="\t"}{gsub("-.*", "", $1); print}' \
   | sort \
  > inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv

# difference between the two
wc -l inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv
wc -l inputOnScratch/gff/${thisspecies}.tsv


# get names of protein from gff file (22,745)
cut -f 1 inputOnScratch/gff/${thisspecies}.tsv > tmp/${thisspecies}-gff-names

# get names of protein from fasta file (26,911)
cut -f 1 -d " " inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv | sort > tmp/${thisspecies}-fasta-names

# get common protein names (22,745)
comm -12 tmp/${thisspecies}-gff-names tmp/${thisspecies}-fasta-names > tmp/${thisspecies}-common-names

# get protein unique to fasta: (4,166)
comm -13 tmp/${thisspecies}-gff-names tmp/${thisspecies}-fasta-names > tmp/${thisspecies}-names-unique-to-fasta


# reduce the file with protein names and length to keep only the common proteins (20,950)
grep -v -f tmp/${thisspecies}-names-unique-to-fasta inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv > inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv

# keep only the proteins that are common to both (20,950)
grep -f tmp/${thisspecies}-common-names inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv > inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered-matching-gff.tsv

# check number of unique protein names: 13898
cut -f1 -d " " inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered-matching-gff.tsv > tmp/protein-names-20950

# keep proteins in gff that are also found in fasta
grep -f tmp/protein-names-20950 inputOnScratch/gff/${thisspecies}.tsv > inputOnScratch/gff/${thisspecies}-filtered.tsv

# this filtering leaves fewer lines: 20950
# wc -l inputOnScratch/gff/${thisspecies}-filtered.tsv
# wc -l inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv

# how much length did this remove? from 10,195,727 to 6,841,174 (a third)
# awk '{SUM+=$2}END{print SUM}' inputOnScratch/sequence-length/${thisspecies}-sequence-length.tsv
# awk '{SUM+=$2}END{print SUM}' inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv

#################################################################################
## Step 3: obtain the longest isoform for each protein
# based on Anindita Brahma's code
module unload compilers
module unload mpi
module load r/recommended

Rscript --vanilla subset-longest-protein-isoform.R inputOnScratch/gff/${thisspecies}-filtered.tsv inputOnScratch/sequence-length/${thisspecies}-sequence-length-filtered.tsv inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt



#################################################################################
## Step 4: subset protein fasta file for just those isoforms
# subset protein sequences for those longest proteins

# update the faa headers for short name (ie just the gene name)
awk '{gsub("-.*", ""); print}' inputOnScratch/${thisspecies}.faa > inputOnScratch/${thisspecies}-short-name.faa

# remove protein sequences that are not the longest isoform
seqtk subseq inputOnScratch/${thisspecies}-short-name.faa \
        inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt \
        > inputOnScratch/primary_transcripts/${thisspecies}-longest-isoforms.faa

# this protein fasta has the longest isoform, with headers containing just gene id (i.e. >XP_014597714.1)

