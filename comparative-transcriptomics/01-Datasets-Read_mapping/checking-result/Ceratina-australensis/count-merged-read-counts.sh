#!/bin/bash

# Calculate merged read counts from FeatureCount output


# Copyright Emeline Favreau, UCL



# files/path needed as input
thispath="/lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/rnaseq-qc/Ceratina-australensis/result/twelve-paired/"

speciespath="/lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/rnaseq-qc/Ceratina-australensis"

# dataset name
datasetname="twelve-paired"

# make a directory for these calculations
mkdir -p tmp/${datasetname}
mkdir -p result/${datasetname}

#sum the merged read count number for each sample
awk 'BEGIN {sum=0} {for (i=3; i<=NF; i++) a[i]+=$i } END {for (i in a) print a[i]}' \
	result/${datasetname}/merged_gene_counts.txt \
	> tmp/${datasetname}/sum-merged_gene_counts

head -n 1 result/${datasetname}/merged_gene_counts.txt \
	| cut -f 3-14 \
	| sed s'/\t/\n/'g \
	| cut -d "_" -f 1 > tmp/${datasetname}/bam-file-names

paste tmp/${datasetname}/bam-file-names \
	tmp/${datasetname}/sum-merged_gene_counts \
	> result/${datasetname}/number-merged-read-counts-in-samples