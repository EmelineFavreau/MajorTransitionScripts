#install.packages("BiocManager")
#BiocManager::install("WGCNA")   # WGCNA is available on CRAN
library(WGCNA)
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator

# load input
load("../05-WGCNA/orthology-dependent/consensusNetwork_orthogroups3718/Consensus-NetworkConstruction-Manual_consensusMEs_MinModuleSize_10.RData")

# follow this
# https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/wgcna.html

allowWGCNAThreads()          # allow multi-threading (optional)

# our input matrix: I cannot find it
input_mat <-
  
  
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)