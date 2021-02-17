#!/bin/bash

# Obtain longest transcript for  each gene in Mellifera
# Copyright Emeline Favreau, UCL

## Analysis overview
### Step 1: get the gene and protein ids from gff files
### Step 2: get the protein length from the proteomes
### Step 3: obtain the longest isoform for each protein
### Step 4: subset protein fasta file for just those isoforms

#################################################################################
## Step 1: get the gene and protein ids from gff files

# obtain gff
mkdir -p inputOnScratch/gff

cd inputOnScratch/gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz

mv GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz Apis_mellifera.gff.gz

thisspecies="Apis_mellifera"

# check that gff is NCBI format (RefSeq)
# zcat ${thisspecies}.gff.gz | head

# this is specific to each GFF (not a norm unless from Ensembl)
# get gene and protein IDs from GFF file, e.g. 118439541 XP_035721480.1
# use uniq to remove duplicated
# useful link : https://www.biostars.org/p/388936/
# remove in column 9 anything that is before and includes ;Dbxref=GeneID:
# remove in column 9 anything that is after ,BEEBASE;
# remove the string ,Genbank:
zcat ${thisspecies}.gff.gz | awk '$3=="CDS"{gsub(".+;Dbxref=GeneID:", "", $9); gsub(",BEEBASE:.+", "", $9); gsub(";Name=.+", "", $9); gsub(",Genbank:", "\t", $9); print $9}' | sort | uniq > ${thisspecies}.tsv

# search for mitochondrial genes
sed -n '/ID=/'p ${thisspecies}.tsv | awk '{ gsub("ID=cds-", ""); gsub(";Parent=.*", ""); print }' > ${thisspecies}-mitochondrial-gene

linestoremove=`wc -l ${thisspecies}-mitochondrial-gene`

# find out the total lines in the tsv file
totallines=`less ${thisspecies}.tsv | wc -l`

# do (b-a)+1 to get the line from where the tsv file should be trimmed
trimfromhere=`expr $totallines - $linestoremove + 1`

# now trim the file
sed -i "${trimfromhere}, ${totallines}d" ${thisspecies}.tsv

cd ..



#################################################################################
## Step 2: get the protein length from the proteomes


mkdir -p sequence-length 

# first, tidy fastaa headers by checking first ( head -n 1 Apis_mellifera.faa | cut -d " " -f 1)
# remove the second part of the sequence header after the gene ID
# fx2tab calculate the sequence length
# sort the result
~/programs/./seqkit fx2tab -n -l ${thisspecies}.faa \
  | awk 'BEGIN{FS="\t"}{gsub(" .*", "", $1); print}' \
  | sort \
  > sequence-length/${thisspecies}-sequence-length.tsv

# remove mitochondrial genes (output non-matching mito file)
grep -v -f gff/Apis_mellifera-mitochondrial-gene sequence-length/${thisspecies}-sequence-length.tsv > sequence-length/${thisspecies}-sequence-length-no-mito.tsv



#################################################################################
## Step 3: obtain the longest isoform for each protein
# based on Anindita's longest_protein_isoform_Step_2.Rmd
module unload compilers
module unload mpi
module load r/recommended

mkdir -p inputOnScratch/longest-isoform

Rscript --vanilla subset-longest-protein-isoform.R inputOnScratch/gff/${thisspecies}.tsv inputOnScratch/sequence-length/${thisspecies}-sequence-length-no-mito.tsv inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt


#################################################################################
## Step 4: subset protein fasta file for just those isoforms
# subset protein sequences for those longest proteins

mkdir -p inputOnScratch/OrthoFinder-input

# update the faa headers for short name (ie just the gene name)
awk '{gsub(" .*", ""); print}' inputOnScratch/${thisspecies}.faa > inputOnScratch/${thisspecies}-short-name.faa

# remove protein sequences that are not the longest isoform
seqtk subseq inputOnScratch/${thisspecies}-short-name.faa \
	inputOnScratch/longest-isoform/${thisspecies}_protein_longest_isoform.txt \
	> inputOnScratch/primary_transcripts/${thisspecies}-longest-isoforms.faa
# this protein fasta has the longest isoform, with headers containing just gene id (i.e. >NP_001010975.1)