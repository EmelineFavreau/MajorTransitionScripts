- Directories minModuleSizexx
  These directories hold the results from the consensus WGCNA meta-analyses performed via
  Z- transform. Files as follows:

--- ModuleTraitRelationships_Heatmap_[species].pdf
	Heatmaps showing the correlations between consensus modules and the traits. Colors
	represent the correlation coefficients with the traits; asterisks denote significance
	of those coefficients at * < 0.05, ** < 0.01, ** < 0.001
	
--- Module-Trait_MetaAnalysis_Results.txt
	Table for all consensus modules showing the meta-analysis Z-score for the traits of
	reproductive / non-reproductive and the associated P-value for that Z-score. Modules 
	with a P < 0.05 were considered significantly associated with traits across the
	species.
	
--- Gene-Module_MetaAnalysis_Results.txt
	Large table of the results of the meta-analysis between genes and modules to identify
	which orthogroups are significantly associated with which consensus modules across
	all species. First two columns for each module color provide the meta-analysis 
	Z-scores and P-values, per orthogroup. Remaining columns per module color provide the 
	individual species' correlation coefficients that went into the meta-analysis.
	Results are only provided for correlations with reproductives because non-reproductive
	correlation coefficients are negative and the associated p-values are identical. This
	table has not been sorted or filtered.
	
--- Gene-Trait_MetaAnalysis_Results.txt
	Large table of the results of the meta-analysis between genes and traits to identify
	which orthogroups are significantly associated with traits across all species. First 
	two columns of the table provide the meta-analysis Z-scores and P-values, per 
	orthogroup. Remaining columns provide the individual species' correlation coefficients
	and their p-values that went into the meta-analysis. Results are only provided for 
	correlations with reproductives because non-reproductive correlation coefficients are 
	negative and associated p-values are identical. This table has not been sorted or 
	filtered.
	
--- orthogroupList_SigModuleAssociated_Only.txt
	List of all orthogroups, filtered for redundancy, identified by meta-analysis to be 
	significantly associated with the consensus modules that are in turn significantly 
	associated with reproductive / non-reproductive trait status. 

--- orthogroupList_SigModule-SigTraitAssociated.txt
	List of all orthogroups, filtered for redundancy, identified by meta-analysis to be 
	significantly associated with the consensus modules that are in turn significantly 
	associated with reproductive / non-reproductive trait status. Further, these 
	orthogroups were also identified by meta-analysis to be significantly associated with
	trait status. This is the most conservative filtering performed. An orthogroup would
	be significantly associated across the majority of species with trait and with trait-
	associated module to be included in this list. 

--- orthogroupList_SigModule-SigTraitAssociated_R.txt
	Orthogroups from the above list (orthogroupList_SigModule-SigTraitAssociated) that
	have a positive meta-analysis Z-score associated with trait, identifying them as 
	significantly correlated with "reproductive" status. List has only been filtered for 
	redundancy.
	
--- orthogroupList_SigModule-SigTraitAssociated_NR.txt
	Orthogroups from the above list (orthogroupList_SigModule-SigTraitAssociated) that
	have a negative meta-analysis Z-score associated with trait, identifying them as 
	significantly correlated with "non-reproductive" status. List has only been filtered 
	for redundancy.
	
	Note: an orthogroup can be included in both an R and an NR list, depending on which 
	module it was correlated with. Some othogroups have meta-analysis Z-scores that are 
	positive for one module and negative for another module. Thus, we only focus on the
	combined list for functional annotation of these genes.

- Individual file: Iterative_MinModuleSize_Results.txt
  Contains a summary of the numbers of modules and genes per module identified by 
  iteratively reducing the minimum module size during manual network construction.
  
- Individual file: consensusNetwork_input_data.Rdata
  An Rdata file of all input files used for manual consensus network construction with 
  WGCNA.
  
- Individual file: Consensus-NetworkConstruction-manual_MinModuleSize_30.RData
  An Rdata file of resulting consensus network using a minimum module size of 30.

- Individual file: Consensus-NetworkConstruction-manual_MinModuleSize_25.RData
  An Rdata file of resulting consensus network using a minimum module size of 25.

- Individual file: Consensus-NetworkConstruction-manual_MinModuleSize_20.RData
  An Rdata file of resulting consensus network using a minimum module size of 20.

- Individual file: Consensus-NetworkConstruction-manual_MinModuleSize_15.RData
  An Rdata file of resulting consensus network using a minimum module size of 15.

- Individual file: Consensus-NetworkConstruction-manual_MinModuleSize_10.RData
  An Rdata file of resulting consensus network using a minimum module size of 10.
  
- Individual file: SampleClustering.pdf
  QC of the clustering of each sample per species in the input to verify that there are
  no outliers.
  


 