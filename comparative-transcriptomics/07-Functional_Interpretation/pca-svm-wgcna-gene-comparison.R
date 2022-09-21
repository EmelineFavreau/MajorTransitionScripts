#---
## Objectives: Answering reviewer comments: 
# which genes are common to SVM, WGCNA and PC9 and PC10?
# author: "E Favreau, UCL"
#date: "2022-09-21"

#---
# aim 1: find the genes that have highest row variance in PC9 and PC10
# aim 2: find common genes between SVM, WGCNA and PC9-PC10
# aim 3: interpret results


###############################################################################
# aim 1: find the genes that have highest row variance
library(tidyverse)
library(tidyr)

# vsdata1 is Variant-stabilized transformed read counts of 
# 3,718 nearly single-copy orthologs across the six species
load("dge-DESeq2/orthology-dependent/orthogroups3718_all6spp_together/VST.Rdata")

# here we look at the top 500 genes (selected by highest row variance)
ntop=500

# calculate the variance for each gene
rv <- rowVars(assay(vsdata1))

# select the ntop genes by variance
select <- order(rv,
                decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsdata1)[select,]))


# find PC9 and PC10 in the matrix of variable loadings 
# (i.e., a matrix whose columns contain the eigenvectors). 
# The function princomp returns this in the element loadings.
# rownames are genes
myfocus <- pca$rotation[ , c("PC9", "PC10")] # OG0006372 -0.173370071  0.042898429

# save absolute values
aload <- abs(myfocus)

# yields the proportional contribution to the each principal component
gene_contribution_pc9pc10 <- sweep(aload, 2, colSums(aload), "/")

# check data
head(gene_contribution_pc9pc10)

# save table for comparison with 71 genetic tool genes and 358 consensus WGCNA genes
write.table(x = gene_contribution_pc9pc10,
            file = "dge-DESeq2/orthology-dependent/orthogroups3718_all6spp_together/orthogroup-proportional-contribution-to-PC9-PC10",
            quote = FALSE,
            sep = "\t")

###############################################################################
# aim 2: find common genes between SVM, WGCNA and PC9-PC10
# list of consensus WGCNA genes (one column, no header, e.g. OG0000394) n = 358
wgcna_genes <- read.table("WGCNA/orthology-dependent/consensusNetwork_orthogroups3718/minModuleSize10/orthogroupList_SigModule-SigTraitAssociated.txt",
                          quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(wgcna_genes) <- "Orthogroup"

# list of 71 genes  (one column, header "Orthogroup", e.g. OG0000706) n = 71
svm_genes <- read.table("svm/result/consensus-wgcna-species-normalised-svm-overlap/overlapping_71_genes.txt",
                          quote="\"", comment.char="", stringsAsFactors=FALSE,
                        header = TRUE)

# any SVM gene in the PC9 or PC10? None
svm_genes$Orthogroup %in% rownames(gene_contribution_pc9pc10)

# any wgcna gene in the PC9 or PC10? one: "OG0002082"
wgcna_genes$Orthogroup[wgcna_genes$Orthogroup %in% rownames(gene_contribution_pc9pc10)]

# checking proportional contribution to each PC
gene_contribution_pc9pc10[rownames(gene_contribution_pc9pc10) == "OG0002082", ]
# PC9        PC10 
# 0.007686782 0.005179726 

# checking these values against the distribution
subset(gene_contribution_pc9pc10, select = "PC9") %>% summary()
# PC9           
# Min.   :7.088e-06  
# 1st Qu.:7.098e-04  
# Median :1.410e-03  
# Mean   :2.000e-03  
# 3rd Qu.:2.663e-03  
# Max.   :1.993e-02  

# which quantile is 0.007686782 in the 9th decile 
# (10% of the data with the highest contribution to PC9)
subset(gene_contribution_pc9pc10, select = "PC9") %>% 
  quantile(probs = seq(.1, .9, by = .1))

subset(gene_contribution_pc9pc10, select = "PC10") %>% summary()
# PC10          
# Min.   :7.968e-06  
# 1st Qu.:6.456e-04  
# Median :1.280e-03  
# Mean   :2.000e-03  
# 3rd Qu.:2.359e-03  
# Max.   :2.990e-02  

# which quantile is 0.005179726  in the 9th decile 
# (10% of the data with the highest contribution to PC10)
subset(gene_contribution_pc9pc10, select = "PC10") %>% 
  quantile(probs = seq(.1, .9, by = .1))

###############################################################################
# aim 3: interpret results

# Only one gene is common between PCA (PC9 and PC10) and consensus WGCNA: "OG0002082"
# Its function is: dopamine N-acetyltransferase-like (XP_033333682.1)
# in Drosophila: "sbr small bristles" (Enables protein N-terminus binding activity)
# Its proportional contributions to PC9 and PC10 are in the top 10% 
# (respectively 0.007686782 and 0.005179726, out of 500 genes with highest variance).
# 


