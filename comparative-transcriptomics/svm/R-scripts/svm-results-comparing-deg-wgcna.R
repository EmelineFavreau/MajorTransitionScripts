library(UpSetR)
library(tidyverse)

# common predictor genes
predictor_genes <- read.csv("result/20-iterations/132_common_predictor_genes",
                            sep = "", stringsAsFactors = FALSE)

# hash tables for each species
Polistes_dominula_hash_table <- read.delim("../orthology-analysis/result/Polistes_dominula_3718_gene_orthogroups_list",
                                           header=FALSE, stringsAsFactors = FALSE)

Polistes_canadensis_hash_table <- read.delim("../orthology-analysis/result/Polistes_canadensis_3718_gene_orthogroups_list",
                                           header=FALSE, stringsAsFactors = FALSE)

Liostenogaster_flavolineata_hash_table <- read.delim("../orthology-analysis/result/Liostenogaster_flavolineata_3718_gene_orthogroups_list",
                                           header=FALSE, stringsAsFactors = FALSE)

Ceratina_australensis_hash_table <- read.delim("../orthology-analysis/result/Ceratina_australensis_3718_gene_orthogroups_list",
                                             header=FALSE, stringsAsFactors = FALSE)

Megalopta_genalis_hash_table <- read.delim("../orthology-analysis/result/Megalopta_genalis_3718_gene_orthogroups_list",
                                                     header=FALSE, stringsAsFactors = FALSE)

Ceratina_calcarata_hash_table <- read.delim("../orthology-analysis/result/Ceratina_calcarata_3718_gene_orthogroups_list",
                                               header=FALSE, stringsAsFactors = FALSE)

# DEGs 
Ceratina_australensis_DEGs <- read.table("../dge-DESeq2/full_set/Ceratina_australensis/Ceratina_australensis_all_DEGs.txt",
                                         quote="\"", comment.char="", stringsAsFactors = FALSE)

Ceratina_calcarata_DEGs <- read.table("../dge-DESeq2/full_set/Ceratina_calcarata/Ceratina_calcarata_all_DEGs.txt",
                                         quote="\"", comment.char="", stringsAsFactors = FALSE)

Megalopta_genalis_DEGs <- read.table("../dge-DESeq2/full_set/Megalopta_genalis/Megalopta_genalis_all_DEGs.txt",
                                         quote="\"", comment.char="", stringsAsFactors = FALSE)

Liostenogaster_flavolineata_DEGs <- read.table("../dge-DESeq2/full_set/Liostenogaster_flavolineata/Liostenogaster_flavolineata_all_DEGs.txt",
                                      quote="\"", comment.char="", stringsAsFactors = FALSE)

Polistes_canadensis_DEGs <- read.table("../dge-DESeq2/full_set/Polistes_canadensis/Polistes_canadensis_all_DEGs.txt",
                                     quote="\"", comment.char="", stringsAsFactors = FALSE)

Polistes_dominula_DEGs <- read.table("../dge-DESeq2/full_set/Polistes_dominula/Polistes_dominula_all_DEGs.txt",
                                               quote="\"", comment.char="", stringsAsFactors = FALSE)

# WGCNA
Polistes_dominula_WGCNAs <- read.table("../WGCNA/full_set/Polistes_dominula/GO_Gene_Lists/Polistes_dominula_WGCNA_Significant_Module_Genes.txt",
                                       quote="\"", comment.char="", stringsAsFactors = FALSE)

Polistes_canadensis_WGCNAs <- read.table("../WGCNA/full_set/Polistes_canadensis/GO_Gene_Lists/Polistes_canadensis_WGCNA_Significant_Module_Genes.txt",
                                       quote="\"", comment.char="", stringsAsFactors = FALSE)

Liostenogaster_flavolineata_WGCNAs <- read.table("../WGCNA/full_set/Liostenogaster_flavolineata/GO_Gene_Lists/Liostenogaster_flavolineata_WGCNA_Significant_Module_Genes.txt",
                                       quote="\"", comment.char="", stringsAsFactors = FALSE)

Ceratina_calcarata_WGCNAs <- read.table("../WGCNA/full_set/Ceratina_calcarata/GO_Gene_Lists/Ceratina_calcarata_WGCNA_Significant_Module_Genes.txt",
                                         quote="\"", comment.char="", stringsAsFactors = FALSE)

Megalopta_genalis_WGCNAs <- read.table("../WGCNA/full_set/Megalopta_genalis/GO_Gene_Lists/Megalopta_genalis_WGCNA_Significant_Module_Genes.txt",
                                                 quote="\"", comment.char="", stringsAsFactors = FALSE)

Ceratina_australensis_WGCNAs <- read.table("../WGCNA/full_set/Ceratina_australensis/GO_Gene_Lists/Ceratina_australensis_WGCNA_Significant_Module_Genes.txt",
                                        quote="\"", comment.char="", stringsAsFactors = FALSE)

##############################################################################
# comparative analysis of orthogroups that are both differentially expressed in each species AND in the common set from SVM
# make a vector for each species with the OG-orthogroups that are DEGs and SVMs
DEGs_SVM_listInput <- list(
  Ceratina_australensis_DEGs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
  Ceratina_australensis_hash_table$V1[Ceratina_australensis_hash_table$V2 %in% Ceratina_australensis_DEGs$V1]],
                           
  Ceratina_calcarata_DEGs_SVM    = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
  Ceratina_calcarata_hash_table$V1[Ceratina_calcarata_hash_table$V2 %in% Ceratina_calcarata_DEGs$V1]],
                           
  Megalopta_genalis_DEGs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
  Megalopta_genalis_hash_table$V1[Megalopta_genalis_hash_table$V2 %in% Megalopta_genalis_DEGs$V1]],
  
  Polistes_canadensis_DEGs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
  Polistes_canadensis_hash_table$V1[Polistes_canadensis_hash_table$V2 %in% Polistes_canadensis_DEGs$V1]],
  
  Liostenogaster_flavolineata_DEGs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
  Liostenogaster_flavolineata_hash_table$V1[Liostenogaster_flavolineata_hash_table$V2 %in% Liostenogaster_flavolineata_DEGs$V1]],
  
  Polistes_dominula_DEGs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
  Polistes_dominula_hash_table$V1[Polistes_dominula_hash_table$V2 %in% Polistes_dominula_DEGs$V1]]
)

# str(DEGs_SVM_listInput)

pdf("DEGs_SVM_comparative_plot.pdf")
# plot comparative analysis of species DEGs in SVM core genes
upset(fromList(DEGs_SVM_listInput),
      sets = c("Ceratina_australensis_DEGs_SVM",
                "Ceratina_calcarata_DEGs_SVM",
               "Megalopta_genalis_DEGs_SVM",
               "Liostenogaster_flavolineata_DEGs_SVM",
               "Polistes_dominula_DEGs_SVM",
               "Polistes_canadensis_DEGs_SVM"),
      order.by = "degree", keep.order = TRUE)

dev.off()
##############################################################################
# comparative analysis of orthogroups that are both in sig WGCNA modules AND in the common set from SVM
# make a vector for each species with the OG-orthogroups that are in sig WGCNA modules and SVMs

WGCNAs_SVM_listInput <- list(
  Ceratina_australensis_WGCNAs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
                                                                    Ceratina_australensis_hash_table$V1[Ceratina_australensis_hash_table$V2 %in% Ceratina_australensis_WGCNAs$V1]],
  
  Ceratina_calcarata_WGCNAs_SVM    = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
                                                                    Ceratina_calcarata_hash_table$V1[Ceratina_calcarata_hash_table$V2 %in% Ceratina_calcarata_WGCNAs$V1]],
  
  Megalopta_genalis_WGCNAs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
                                                                Megalopta_genalis_hash_table$V1[Megalopta_genalis_hash_table$V2 %in% Megalopta_genalis_WGCNAs$V1]],
  
  Polistes_canadensis_WGCNAs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
                                                                  Polistes_canadensis_hash_table$V1[Polistes_canadensis_hash_table$V2 %in% Polistes_canadensis_WGCNAs$V1]],
  
  Liostenogaster_flavolineata_WGCNAs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
                                                                          Liostenogaster_flavolineata_hash_table$V1[Liostenogaster_flavolineata_hash_table$V2 %in% Liostenogaster_flavolineata_WGCNAs$V1]],
  
  Polistes_dominula_WGCNAs_SVM = predictor_genes$predictor_gene[predictor_genes$predictor_gene %in% 
                                                                Polistes_dominula_hash_table$V1[Polistes_dominula_hash_table$V2 %in% Polistes_dominula_WGCNAs$V1]]
)

# str(WGCNAs_SVM_listInput)
pdf("WGCNAs_SVM_comparative_plot.pdf")
# plot comparative analysis of species WGCNAs in SVM core genes
upset(fromList(WGCNAs_SVM_listInput),
      sets = c("Ceratina_australensis_WGCNAs_SVM",
               "Ceratina_calcarata_WGCNAs_SVM",
               "Megalopta_genalis_WGCNAs_SVM",
               "Liostenogaster_flavolineata_WGCNAs_SVM",
               "Polistes_dominula_WGCNAs_SVM",
               "Polistes_canadensis_WGCNAs_SVM"),
      order.by = "degree", keep.order = TRUE)
dev.off()
