# running DESeq2 on data (using resampling approach)
scripts/desq2_fulldata_2021.03.06.Rmd

# making a heatmap of Variance Stabilized read counts
scripts/2021-08-03-heatmap-variant-stabilised-read-count.R

# Orthology-dependent (assumes shared gene history)
## 1. code
### 1.1 deseq2_orthogroups_all6spp+bees+wasps.Rmd
Performs the ortholog-dependent combined DESeq2 method as in Morandin et al. 2017, including VST [combined or separate]. Unlike with full gene set, this program combines the data for all species in the DEG analysis using species and clade as additional variables. No resampling performed here. Results from these analyses were reported in the paper. 

### 1.2 deseq2_orthogroups3718.Rmd
Performs an ortholog-dependent DESeq2 analysis to produce an ortholog-dependent VST for each species for WGCNA. It will also do N = 3 resampling if needed. This code was only used to generate VST for WGCNA.


## 2. results
### 2.1 Species' DEG results directory
**2.1.2 [species]_all6spp_all_DEGs.txt**
	Contains the list of all differentially expressed genes (DEGs) identified for a given species using the cutoff parameter of FDR < 0.05 (after restricting to orthologs).
**[species]_all6spp_DESeq2_results.txt**
	The log-fold expression, p-value, and FDR-adjusted p-value for all genes in the analysis (after restricting to orthologs).
**[species]_all6spp_upNR_DEGs.txt**
	The list of all DEGs identified as upregulated in non-reproductive phenotype (after restricting to orthologs).
**[species]_all6spp_upR_DEGs.txt**
	The list of all DEGs identified as upregulated in reproductive phenotype (after restricting to orthologs).


# Species-specific (orthology-agnostic)
## 1. code
### 1.1 deseq2_fulldata_2021.05.28.Rmd**
This file will perform the ortholog-agnostic DESeq2 analysis, including variance stabilized transform (VST) and N = 3 resampling if needed. The resampling is a switch that needs to be set at the top of the file.  Performs analyses on all species simultaneously.

## 1.2 deseq_resample_counter.pl
Counts the number of times a gene is flagged as significant by DESeq2 iteratively in order to calculate how often resampling yields the same genes. Must be run in order to complete the resampling analysis and before 1.3.		

## 1.3 DEG_ComparePointEst2Resampling.Rmd
Takes the results from the perl script in 1.2 and compares it to the point estimated DEGs identified from 1.1. Generates p-values reported in the paper.

## 2. results
### 2.1 Species' DEG results directory
**2.1.2 [species]_all_DEGs.txt**
	Contains the list of all differentially expressed genes (DEGs) identified for a given species using the cutoff parameter of FDR < 0.05.
**[species]_DESeq2_results.txt**
	The log-fold expression, p-value, and FDR-adjusted p-value for all genes in the analysis. 
**[species]_upNR_DEGs.txt**
	The list of all DEGs identified as upregulated in non-reproductive phenotype.
**[species]_upR_DEGs.txt**
	The list of all DEGs identified as upregulated in reproductive phenotype.
### 2.2 DEGs_Summary_All_Species.txt
Summary statistics of DEGs identified for all species included in the analysis. 
### 2.3 Raw_Read_Counts_Summary_Statistics.txt
Summary statistics of raw read counts used for the DEG analysis for all species. 

