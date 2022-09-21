#---
## Objectives: Answering reviewer 2 comment 2: 
# Draw distribution curve of PC9 and PC10 loading
# author: "E Favreau, UCL"
#date: "2022-09-21"

#---
# aim 1: shuffle the reproductive/non-reproductive labels for the data many times 
# aim 2: report the distribution of data when shuffled
# aim 3: interpret results


###############################################################################
# aim 1: shuffle the reproductive/non-reproductive labels for the data many times 
library(tidyverse)
library(tidyr)

# list of 500 genes that are top contributor to PC9 and PC10
gene_contribution_pc9pc10 <-  read.table("dge-DESeq2/orthology-dependent/orthogroups3718_all6spp_together/orthogroup-proportional-contribution-to-PC9-PC10",
                                         quote="\"", comment.char="", stringsAsFactors=FALSE)

# subset to the top 10 genes: the focus of our simulation
top10_gene_contributing_pc9 <- rownames(gene_contribution_pc9pc10)[
  order(gene_contribution_pc9pc10$PC9, decreasing = T)][1:10]

# proportion of each of them
top10_gene_contribution_pc9 <- 
  gene_contribution_pc9pc10[rownames(gene_contribution_pc9pc10) %in% 
                              top10_gene_contributing_pc9,]

# vsdata1 is Variant-stabilized transformed read counts of 
# 3,718 nearly single-copy orthologs across the six species
load("dge-DESeq2/orthology-dependent/orthogroups3718_all6spp_together/VST.Rdata")

# here we look at the top 500 genes (selected by highest row variance)
ntop=500

# run the randomisation 10 times
random_times <- 10


# instantiate a list for result
randomlist <- list()

# for each iteration
for(i in 1:random_times){
  # copy data for this randomisation
  random_vsdata <- vsdata1 
  
  # randomise the phenotype
  random_vsdata$R.NR_Phenotype <- sample(random_vsdata$R.NR_Phenotype)
  
  # calculate the variance for each gene
  rv <- rowVars(assay(random_vsdata))
  
  # select the ntop genes by variance
  select <- order(rv,
                  decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(random_vsdata)[select,]))
  
  # what are the loadings of PC9 and PC10?
  # find PC9 and PC10 in the matrix of variable loadings 
  # rownames are genes
  myfocus <- pca$rotation[ , c("PC9", "PC10")] 
  
  # save absolute values
  aload <- abs(myfocus)
  
  # yields the proportional contribution to the each principal component
  randomgene_contribution_pc9pc10 <- sweep(aload, 2, colSums(aload), "/")
  
  # find our 10 focus genes
  randomdata <- randomgene_contribution_pc9pc10[rownames(randomgene_contribution_pc9pc10) %in% 
                                    top10_gene_contributing_pc9,]
  
  # store it in a list
  randomlist[[i]] <- randomdata
}

###############################################################################
# aim 2: report the distribution of data when  phenotype is shuffled

# tidy the list into a data ready for plot
randomlist[[1]][1:10] #PC9
rownames(randomlist[[1]])
#randomlist[[1]][11:20] #PC10
data <- data.frame(
orthogroup = c(rownames(randomlist[[1]]),
         rownames(randomlist[[2]]),
         rownames(randomlist[[3]]),
         rownames(randomlist[[4]]),
         rownames(randomlist[[5]]),
         rownames(randomlist[[6]]),
         rownames(randomlist[[7]]),
         rownames(randomlist[[8]]),
         rownames(randomlist[[9]]),
         rownames(randomlist[[10]])),

proportional_contribution_to_PC9 = c(randomlist[[1]][1:10],
          randomlist[[2]][1:10],
          randomlist[[3]][1:10],
          randomlist[[4]][1:10],
          randomlist[[5]][1:10],
          randomlist[[6]][1:10],
          randomlist[[7]][1:10],
          randomlist[[8]][1:10],
          randomlist[[9]][1:10],
          randomlist[[10]][1:10]))


# Plot PC9
data %>%
  # randomise data in black
  ggplot( aes(x=orthogroup, y=proportional_contribution_to_PC9, fill=orthogroup)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  
  # real data in red
  geom_point(data = data.frame(orthogroup = rownames(top10_gene_contribution_pc9),
                               proportional_contribution_to_PC9 = top10_gene_contribution_pc9$PC9),
             aes(x=orthogroup, y=proportional_contribution_to_PC9),
             color = 'red') +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 45, hjust=1)) +
  ggtitle("Black Loadings: Random Phenotype; Red Loadings: Real Phenotype") +
  xlab("Orthogroup")

# save figure for suppl
ggsave("../../MajorTransitionScripts/comparative-transcriptomics/Figures/Supplementary/2022-09-21-PCA-randomisation.tiff")


###############################################################################
# aim 3: interpret results