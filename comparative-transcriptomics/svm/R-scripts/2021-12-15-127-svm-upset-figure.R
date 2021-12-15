# This script makes the upset plot for the manuscript.
# get libraries
basic_libraries <- c("UpSetR",
                     "tidyverse")

for (lib in basic_libraries) {
  if (require(package = lib, character.only = TRUE)) {
    print("Successful")
  } else {
    print("Installing")
    install.packages(lib)
    library(lib, character.only = TRUE )
  }
}

# common predictor genes (colname predictor_gene)
predictor_genes <- 
  read.csv("../result/species-normalised/regression-on-all/127_common_predictor_genes",
                            sep = "", stringsAsFactors = FALSE)

# list of all orthogroups (colname V1, V2)
all_orthogroups_list <- 
  read.delim("../input/Ceratina_australensis_3718_gene_orthogroups_list",
                             header=FALSE, stringsAsFactors=FALSE)

# hash tables for each species (colname V1, V2)
Polistes_dominula_hash_table <- 
  read.delim("../../orthology-analysis/result/Polistes_dominula_3718_gene_orthogroups_list",
   header=FALSE, stringsAsFactors = FALSE)

Polistes_canadensis_hash_table <- 
  read.delim("../../orthology-analysis/result/Polistes_canadensis_3718_gene_orthogroups_list",
  header=FALSE, stringsAsFactors = FALSE)

Liostenogaster_flavolineata_hash_table <- 
  read.delim("../../orthology-analysis/result/Liostenogaster_flavolineata_3718_gene_orthogroups_list",
   header=FALSE, stringsAsFactors = FALSE)

Ceratina_australensis_hash_table <- 
  read.delim("../../orthology-analysis/result/Ceratina_australensis_3718_gene_orthogroups_list",
  header=FALSE, stringsAsFactors = FALSE)

Megalopta_genalis_hash_table <- 
  read.delim("../../orthology-analysis/result/Megalopta_genalis_3718_gene_orthogroups_list",
    header=FALSE, stringsAsFactors = FALSE)

Ceratina_calcarata_hash_table <- 
  read.delim("../../orthology-analysis/result/Ceratina_calcarata_3718_gene_orthogroups_list",
     header=FALSE, stringsAsFactors = FALSE)

# list of predictor genes per species (colname V1)
Ceratina_australensis_predicted_gene_vec <- 
  read.table("../result/species-normalised/regression-on-all/Ceratina_australensis_svm_predicted_gene_list_tidy",
    quote="\"", comment.char="", stringsAsFactors = FALSE)

Ceratina_calcarata_predicted_gene_vec <- 
  read.table("../result/species-normalised/regression-on-all/Ceratina_calcarata_svm_predicted_gene_list_tidy",
    quote="\"", comment.char="", stringsAsFactors = FALSE)

Megalopta_genalis_predicted_gene_vec <- 
  read.table("../result/species-normalised/regression-on-all/Megalopta_genalis_svm_predicted_gene_list_tidy",
     quote="\"", comment.char="", stringsAsFactors = FALSE)

Polistes_canadensis_predicted_gene_vec <- 
  read.table("../result/species-normalised/regression-on-all/Polistes_canadensis_svm_predicted_gene_list_tidy",
      quote="\"", comment.char="", stringsAsFactors = FALSE)

Polistes_dominula_predicted_gene_vec <- 
  read.table("../result/species-normalised/regression-on-all/Polistes_dominula_svm_predicted_gene_list_tidy",
      quote="\"", comment.char="", stringsAsFactors = FALSE)

Liostenogaster_flavolineata_predicted_gene_vec <- 
  read.table("../result/species-normalised/regression-on-all/Liostenogaster_flavolineata_svm_predicted_gene_list_tidy",
     quote="\"", comment.char="", stringsAsFactors = FALSE)


# Aim 1: Overlap of SVM predictor genes

# all_orthogroups_list$V1[all_orthogroups_list$V1 %in% Ceratina_australensis_predicted_gene_vec$V1]

# comparative analysis of orthogroups that are predictro genes in each svm
SVM_listInput <- list(
  Ceratina_australensis_SVM       = all_orthogroups_list$V1[all_orthogroups_list$V1 %in% Ceratina_australensis_predicted_gene_vec$V1],
  Ceratina_calcarata_SVM          = all_orthogroups_list$V1[all_orthogroups_list$V1 %in% Ceratina_calcarata_predicted_gene_vec$V1],
  Megalopta_genalis_SVM           = all_orthogroups_list$V1[all_orthogroups_list$V1 %in% Megalopta_genalis_predicted_gene_vec$V1],
  Polistes_canadensis_SVM         = all_orthogroups_list$V1[all_orthogroups_list$V1 %in% Polistes_canadensis_predicted_gene_vec$V1],
  Polistes_dominula_SVM           = all_orthogroups_list$V1[all_orthogroups_list$V1 %in% Polistes_dominula_predicted_gene_vec$V1],
  Liostenogaster_flavolineata_SVM = all_orthogroups_list$V1[all_orthogroups_list$V1 %in% Liostenogaster_flavolineata_predicted_gene_vec$V1])

# str(SVM_listInput)

UpSetR::upset(fromList(SVM_listInput),
      sets = c("Ceratina_australensis_SVM",
               "Ceratina_calcarata_SVM",
               "Megalopta_genalis_SVM",
               "Liostenogaster_flavolineata_SVM",
               "Polistes_dominula_SVM",
               "Polistes_canadensis_SVM"),
      order.by = "degree", keep.order = TRUE)

upset(fromList(SVM_listInput),
      sets = c("Ceratina_australensis_SVM",
               "Ceratina_calcarata_SVM",
               "Megalopta_genalis_SVM",
               "Liostenogaster_flavolineata_SVM",
               "Polistes_dominula_SVM",
               "Polistes_canadensis_SVM"),
      order.by = "freq")

# keep a copy as pdf
pdf("SVM_3718_Upset_comparative_plot_species_normalised.pdf")
upset(fromList(SVM_listInput),
      sets = c("Ceratina_australensis_SVM",
               "Ceratina_calcarata_SVM",
               "Megalopta_genalis_SVM",
               "Liostenogaster_flavolineata_SVM",
               "Polistes_dominula_SVM",
               "Polistes_canadensis_SVM"),
      order.by = "degree", keep.order = TRUE)
dev.off()

#install.packages('ComplexUpset')
library("ComplexUpset")

# upset(fromList(SVM_listInput),
#       sets = c("Ceratina_australensis_SVM",
#                "Ceratina_calcarata_SVM",
#                "Megalopta_genalis_SVM",
#                "Liostenogaster_flavolineata_SVM",
#                "Polistes_dominula_SVM",
#                "Polistes_canadensis_SVM"),
#       order.by = "degree",
#       keep.order = TRUE,
#       queries=list(
#         upset_query(
#           intersect=c('Ceratina_australensis_SVM',
#                       'Ceratina_calcarata_SVM',
#                       'Megalopta_genalis_SVM'),
#           color='red',
#           fill='red',
#           only_components=c('intersections_matrix', 'Intersection size')
#         )))


### learning new package
library(ggplot2)
library(ComplexUpset)
library(ggbeeswarm)

# if(!require(ggplot2movies)) install.packages('ggplot2movies')
# movies = ggplot2movies::movies
# genres = c('Action', 'Animation', 'Comedy', 'Drama', 'Documentary', 'Romance')
# 
# upset(
#   movies,
#   genres,
#   
#   queries=list(
#     upset_query(
#       intersect=c('Drama', 'Comedy'),
#       color='red',
#       fill='red',
#       only_components=c('intersections_matrix', 'Intersection size')
#     ),
#     upset_query(
#       set='Drama',
#       fill='blue'
#     )
#   ),
#   min_size=10,
#   width_ratio=0.1
# )

# testing with my data
orthogroups = SVM_listInput
species = c("Ceratina_australensis",
            "Ceratina_calcarata",
            "Megalopta_genalis",
            "Liostenogaster_flavolineata",
            "Polistes_dominula",
            "Polistes_canadensis")


# need a dataframe with binary columns representing membership in classes
# orthogroups in row.names

colouredUpsetdf <- data.frame(row.names = all_orthogroups_list$V1)
# species in columns
colouredUpsetdf$Ceratina_australensis       <- 0
colouredUpsetdf$Ceratina_calcarata          <- 0
colouredUpsetdf$Polistes_dominula           <- 0
colouredUpsetdf$Polistes_canadensis         <- 0
colouredUpsetdf$Megalopta_genalis           <- 0
colouredUpsetdf$Liostenogaster_flavolineata <- 0

# add 1 if the orthogroup is present in the species
for(i in all_orthogroups_list$V1){
  if(i %in% orthogroups$Ceratina_australensis_SVM){
    colouredUpsetdf$Ceratina_australensis[row.names(colouredUpsetdf) == i] <- 1
  }
  
  if(i %in% orthogroups$Ceratina_calcarata_SVM){
    colouredUpsetdf$Ceratina_calcarata[row.names(colouredUpsetdf) == i] <- 1
  }
  
  if(i %in% orthogroups$Polistes_dominula_SVM){
    colouredUpsetdf$Polistes_dominula[row.names(colouredUpsetdf) == i] <- 1
  }
  
  if(i %in% orthogroups$Polistes_canadensis_SVM){
    colouredUpsetdf$Polistes_canadensis[row.names(colouredUpsetdf) == i] <- 1
  }
  
  if(i %in% orthogroups$Megalopta_genalis_SVM){
    colouredUpsetdf$Megalopta_genalis[row.names(colouredUpsetdf) == i] <- 1
  }
  
  if(i %in% orthogroups$Liostenogaster_flavolineata_SVM){
    colouredUpsetdf$Liostenogaster_flavolineata[row.names(colouredUpsetdf) == i] <- 1
  }
  
}
  
# check each species count
colSums(colouredUpsetdf)

# a vector of blues for the bees
blues_vec <- c("#9ecae1",
               "#6baed6",
               "#3182bd",
               "#08519c",
               "#deebf7")

# a vector of oranges for the wasps
oranges_vec <- c("#fdae6b",
                 "#fd8d3c",
                 "#f16913",
                 "#d94801","#fff5eb")

# a table with taxon level for each species
taxa_metadata <- data.frame(
  set = factor(species, levels = species),
  taxatype = c("bee", "bee","wasp",         
            "bee", "wasp", "wasp")
)

# plot and save a pdf version of the plot
pdf("SVM_3718_ComplexUpset_comparative_plot_species_normalised.pdf")
ComplexUpset::upset(
  colouredUpsetdf,
  species,
  queries=list(
    upset_query(
      intersect=c("Polistes_canadensis"),
      color=oranges_vec[1],
      fill=oranges_vec[1],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Polistes_dominula"),
      color=oranges_vec[2],
      fill=oranges_vec[2],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Liostenogaster_flavolineata"),
      color=oranges_vec[3],
      fill=oranges_vec[3],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Megalopta_genalis"),
      color=blues_vec[1],
      fill=blues_vec[1],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Ceratina_calcarata"),
      color=blues_vec[2],
      fill=blues_vec[2],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Ceratina_australensis"),
      color=blues_vec[3],
      fill=blues_vec[3],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Polistes_canadensis",
                  "Polistes_dominula",
                  "Liostenogaster_flavolineata"),
      color=oranges_vec[4],
      fill=oranges_vec[4],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
    intersect=c("Polistes_canadensis",
                  "Polistes_dominula"),
      color=oranges_vec[4],
      fill=oranges_vec[4],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Polistes_canadensis",
                  "Liostenogaster_flavolineata"),
      color=oranges_vec[4],
      fill=oranges_vec[4],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Polistes_dominula",
                  "Liostenogaster_flavolineata"),
      color=oranges_vec[4],
      fill=oranges_vec[4],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Ceratina_australensis",
                  "Ceratina_calcarata",
                  "Megalopta_genalis"),
      color=blues_vec[4],
      fill=blues_vec[4],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Ceratina_australensis",
                  "Ceratina_calcarata"),
      color=blues_vec[4],
      fill=blues_vec[4],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Ceratina_australensis",
                  "Megalopta_genalis"),
      color=blues_vec[4],
      fill=blues_vec[4],
      only_components=c('intersections_matrix', 'Intersection size')),
    upset_query(
      intersect=c("Ceratina_calcarata",
                  "Megalopta_genalis"),
      color=blues_vec[4],
      fill=blues_vec[4],
      only_components=c('intersections_matrix', 'Intersection size'))) ,
  max_size=300,
  width_ratio=0.1,sort_intersections=FALSE,
  intersections=list(
    c("Ceratina_australensis",
      "Ceratina_calcarata",
      "Megalopta_genalis",
      "Liostenogaster_flavolineata",
      "Polistes_dominula",
      "Polistes_canadensis"),
    
    c("Ceratina_australensis",
      "Ceratina_calcarata",
      "Megalopta_genalis",
      "Polistes_dominula",
      "Polistes_canadensis"),
    
    c("Ceratina_australensis",
      "Ceratina_calcarata",
      "Megalopta_genalis",
      "Liostenogaster_flavolineata",
      "Polistes_canadensis"),
    
    c("Ceratina_australensis",
      "Ceratina_calcarata",
      "Liostenogaster_flavolineata",
      "Polistes_dominula",
      "Polistes_canadensis"),
    c("Ceratina_australensis",
      "Ceratina_calcarata",
      "Megalopta_genalis",
      "Liostenogaster_flavolineata",
      "Polistes_dominula"),
      
    c("Ceratina_australensis",
      "Megalopta_genalis",
      "Liostenogaster_flavolineata",
      "Polistes_dominula",
      "Polistes_canadensis"),
    c("Ceratina_calcarata",
      "Megalopta_genalis",
      "Liostenogaster_flavolineata",
      "Polistes_dominula",
      "Polistes_canadensis"),



    
    c("Polistes_canadensis",
      "Polistes_dominula",
      "Liostenogaster_flavolineata"),
    
    c("Polistes_canadensis",
      "Polistes_dominula"),
    
    c("Polistes_canadensis",
      "Liostenogaster_flavolineata"),
    
    c("Polistes_dominula",
      "Liostenogaster_flavolineata"),
    c("Ceratina_australensis",
      "Ceratina_calcarata",
      "Megalopta_genalis"),
    c("Ceratina_australensis",
      "Ceratina_calcarata"),
    c("Ceratina_australensis",
      "Megalopta_genalis"),
    c("Ceratina_calcarata",
      "Megalopta_genalis"),
    "Ceratina_australensis",
    "Ceratina_calcarata",
    "Megalopta_genalis",
    "Liostenogaster_flavolineata",
    "Polistes_dominula",
    "Polistes_canadensis"
  ), sort_sets=FALSE,
  
  stripes=upset_stripes(data=taxa_metadata,
    mapping=aes(color=taxa_metadata$taxatype),
    colors=c("bee" = blues_vec[5],
             "wasp" = oranges_vec[5]
             ) 
  )
)
dev.off()

# checking the data are same as plotted
# colouredUpsetdf %>%  filter(total == 1) %>% filter(Polistes_dominula == 1) %>% nrow()