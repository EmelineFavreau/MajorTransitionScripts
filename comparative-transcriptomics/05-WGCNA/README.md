# Consensus

## 1. code
### 1.1 WGCNA_orthogroups_all6spp+bees+wasps_Network_Construction.Rmd
Performs manual construction of the network using the data for bees+wasps, bees only, or wasps only. A switch is set to indicate which analysis is to be performed. Update your directories accordingly. There is an option to perform reshuffling of the trait labels and generate those data as well.

### 1.2 WGCNA_orthogroups_all6spp+bees+wasps_Trait_Correlation.Rmd
Performs the trait correlation analysis for bees+wasps, bees only, or wasps only using the Topological Overlap Matrix and consensus modules identified in 1.1. A switch is set to indicate which analysis is to be performed. Update your directories accordingly.

### 1.3 WGCNA_orthogroups_all6spp+bees+wasps_Trait_Correlation_resampling.Rmd
Performs the trait correlation analysis for bees+wasps, bees only, or wasps only using the Topological Overlap Matrix and consensus modules identified in 1.1 **for the reshuffled data only**. A switch is set to indicate which analysis is to be performed. Update your directories accordingly.

### 1.4 WGCNA_orthogroups_all6spp+bees+wasps_MainFigure.Rmd
Produces the main figure for the paper using the results of the consensus WGCNA, SVM, and DE analyses.

## 2. results
### 2.1 Directories: bees+wasps, bees_only, wasps_only
	Each contains the results for three separate analyses, the bees and wasps combined (N 
	= 3718 orthologs), among bees only (N = 5787) and among wasps only (N = 6983)
**2.1.1 Iterative_MinModuleSize_Results.txt**
	Contains a summary of the numbers of modules and genes per module identified by 		
	iteratively reducing the minimum module size during manual network construction.
**2.1.2 SampleClustering.pdf**
    QC of the clustering of each sample per species in the input to verify that there are 
    no outliers.
**2.1.3 Directories: minModuleSize[xx]**
    These directories hold the results from the consensus WGCNA meta-analyses performed 
    via Z-transform. xx refers to the minimum module size
**2.1.3.1 ModuleTraitRelationships_Heatmap_[species].pdf**
	Heatmaps showing the correlations between consensus modules and the traits. Colors 
	represent the correlation coefficients with the traits; asterisks denote significance 
	of those coefficients at * < 0.05, ** < 0.01, ** < 0.001
**2.1.3.2 Module-Trait_MetaAnalysis_Results.txt**
	Table for all consensus modules showing the meta-analysis Z-score for the traits of 
	reproductive / non-reproductive and the associated P-value for that Z-score. Modules 
	with a P < 0.05 were considered significantly associated with traits across the 
	species.
**2.1.3.4 Gene-Module_MetaAnalysis_Results.txt**
	Large table of the results of the meta-analysis between genes and modules to identify 
	which orthogroups are significantly associated with which consensus modules across all 
	species. First two columns for each module color provide the meta-analysis Z-scores 
	and P-values, per orthogroup. Remaining columns per module color provide the 
	individual species' correlation coefficients that went into the meta-analysis. Results 
	are only provided for correlations with reproductives because non-reproductive 
	correlation coefficients are negative and the associated p-values are identical. This 
	table has not been sorted or filtered.
**2.1.3.5 Gene-Trait_MetaAnalysis_Results.txt**
	Large table of the results of the meta-analysis between genes and traits to identify w
	hich orthogroups are significantly associated with traits across all species. First 
	two columns of the table provide the meta-analysis Z-scores and P-values, per 
	orthogroup. Remaining columns provide the individual species' correlation coefficients 
	and their p-values that went into the meta-analysis. Results are only provided for 
	correlations with reproductives because non-reproductive correlation coefficients are 
	negative and associated p-values are identical. This table has not been sorted or f
	iltered.
**2.1.3.6 orthogroupList_SigModule-SigTraitAssociated.txt**
	List of all orthogroups, filtered for redundancy, identified by meta-analysis to be 
	significantly associated with the consensus modules that are in turn significantly 
	associated with reproductive / non-reproductive trait status. Further, these 
	orthogroups were also identified by meta-analysis to be significantly associated with 
	trait status. This is the most conservative filtering performed. An orthogroup would 
	be significantly associated across the majority of species with trait and with 
	trait-associated module to be included in this list. 

# Species-specific

## 1. code
### 1.1 WGCNA_species-specific_orthogroup_agnostic
Performs automatic construction of each species' coexpression network using the variance-stabilized transformed values from DESeq2. This analysis assumes nothing about shared gene history of the species, and uses all genes expressed by a given species to construst the network. The result is a gene list of significantly trait-associated genes for each species, which could then be compared for shared functions (with gene ontology analysis) without relying on shared gene ancestry.

## 2. results
Each species' folder contains the following results:

**2.1 [species]_Genes_Module_Correlations.txt**
	Table with the genes and their correlation values to each of the modules identified in 
	the network.
**2.2 [species]_Choose_Soft_Threshold_power.pdf**
	Graph showing the results of choosing the soft thresholding power for that species.
**2.3 [species]_Network_Dendrogram_Signed_Network.pdf**
	Dendrogram of the signed network for that species.
**2.4 [species]_soft_threshold_results.txt**
	Tabular results of choosing the soft thresholding power.
**2.5 [species]_Trait_Correlations_Signed_Network.pdf**
	Heatmap of the trait-associated modules and their correlation values and significance 
	for the network for that species.
**2.6 directory: GO_Gene_Lists**
	The gene lists used for Gene Ontology (GO) analysis. 
**2.6.1 [species]_WGCNA_Significant_Module_Genes.txt**
	All genes expressed in the significant trait-associated modules of that species' 
	network.
**2.6.2 [species]_WGCNA_Sig_Module_DEG_Overlap.txt**
	The genes that overlap between the significantly trait-associated modules and those 
	identified independently by DEG analysis with DESeq2.
**2.6.3 [species]_WGCNA_Input_Genes.txt**
	The genes that went into WGCNA, used as the universe for GO analysis.
	
	
#  scripts
# analysing results with TopGO
scripts/2021-10-11-topgo-on-consensus-wgcna.Rmd
scripts/2021-10-25-go-terms-analysis-consensus-wgcna.Rmd
