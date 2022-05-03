# Heatmap of DEG genes for all species
# DEG = differenrially expressed genes

# import data 
# columns: sample names
# row.names: gene names (need to be translated in orthogroups)
# rows: for each orthogroup, the sample has one variant-stabilised read count
# Caust.v2_010925
Ceratina_australensis_VST <- read.delim("input/Ceratina_australensis_VST.txt",
                                        stringsAsFactors = FALSE)     

# gene-LOC108630428
Ceratina_calcarata_VST <- read.delim("input/Ceratina_calcarata_VST.txt",
                                        stringsAsFactors = FALSE)   

# Lflavo2a009208
Liostenogaster_flavolineata_VST <- read.delim("input/Liostenogaster_flavolineata_VST.txt",
                                        stringsAsFactors = FALSE)   

# gene-LOC117229605
Megalopta_genalis_VST <- read.delim("input/Megalopta_genalis_VST.txt",
                                        stringsAsFactors = FALSE)   

# gene-LOC107072336
Polistes_dominula_VST <- read.delim("input/Polistes_dominula_VST.txt",
                                        stringsAsFactors = FALSE)   
# gene-LOC106793534
Polistes_canadensis_VST <- read.delim("input/Polistes_canadensis_VST.txt",
                                        stringsAsFactors = FALSE)   

# columns
# rows
# 

# exploring results
options(tidyverse.quiet = TRUE)
basic_libraries <- c("ggplot2",
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

# combine info
all_species_table <- rbind(Pdom_table[, 1:2],
                           Pcan_table[, 1:2],
                           Caus_table[, 1:2],
                           Ccal_table[, 1:2],
                           Mgen_table[, 1:2],
                           Lfla_table[, 1:2])

# remove NA
all_species_table <- all_species_table[complete.cases(all_species_table), ]

# add species column
all_species_table$species <- c(rep("Pdom", time =3618),
                               rep("Pcan", time =3618),
                               rep("Caus", time =3618),
                               rep("Ccal", time =3618),
                               rep("Mgen", time =3618),
                               rep("Lfla", time =3618))

#all_species_table$species <- factor(all_species_table$species)

# need to sort the orthogroup by one species?
# all_species_table$orthogroup <- factor(all_species_table$orthogroup,
#  levels=unique(all_species_table$orthogroup[all_species_table$species == "Ccal"]))


# Plot OG in x axis, steps in y axis
# one line per species

# # Plot all species
# all_species_table %>% 
# ggplot( aes(x=orthogroup,
#             y=nfeatures,
#             group = species,
#             colour = species)) +
#   geom_line() 
# 
# # Plot all species, only 1000 bst orthogroups
# all_species_table %>% 
#   filter(species %in% c("Caus", "Ccal", "Mgen")) %>% 
#   ggplot( aes(x=orthogroup,
#               y=nfeatures,
#               group = species,
#               colour = species)) +
#   geom_line() 
# 
# # Plot all species, only 1000 bst orthogroups
# all_species_table %>% 
#   filter(species %in% c("Pdom", "Pcan", "Lfla")) %>% 
#   ggplot( aes(x=orthogroup,
#               y=nfeatures,
#               group = species,
#               colour = species)) +
#   geom_line() 

# number of features kept
Caus <- 554
Ccal <- 621
Mgen <- 548
Pdom <- 619
Pcan <- 594
Lfla <- 516

# subset of genes that are predictor in each model
Caus_genes <- all_species_table %>%
  filter(species == "Caus" & nfeatures <= Caus) %>% select(orthogroup)

Ccal_genes <- all_species_table %>%
  filter(species == "Ccal" & nfeatures <= Ccal) %>% select(orthogroup)

Mgen_genes <- all_species_table %>%
  filter(species == "Mgen" & nfeatures <= Mgen) %>% select(orthogroup)

Pdom_genes <- all_species_table %>%
  filter(species == "Pdom" & nfeatures <= Pdom) %>% select(orthogroup)

Pcan_genes <- all_species_table %>%
  filter(species == "Pcan" & nfeatures <= Pcan) %>% select(orthogroup)

Lfla_genes <- all_species_table %>%
  filter(species == "Lfla" & nfeatures <= Lfla) %>% select(orthogroup)

# list of all orthogroups
orthogroup_vec <- readcounts_3718orthogroups_6species$orthogroup

# predictor genes: from the orignal 3718, remove those that are prior to the threshold
# add those that were never in the output (1:100)


Caus_predictor_genes <- c(as.character(Caus_genes$orthogroup),
                          orthogroup_vec[!orthogroup_vec %in% as.character(Caus_table$orthogroup)])

Ccal_predictor_genes <- c(as.character(Ccal_genes$orthogroup),
                          orthogroup_vec[!orthogroup_vec %in% as.character(Ccal_table$orthogroup)])

Pdom_predictor_genes <- c(as.character(Pdom_genes$orthogroup),
                          orthogroup_vec[!orthogroup_vec %in% as.character(Pdom_table$orthogroup)])

Pcan_predictor_genes <- c(as.character(Pcan_genes$orthogroup),
                          orthogroup_vec[!orthogroup_vec %in% as.character(Pcan_table$orthogroup)])

Mgen_predictor_genes <- c(as.character(Mgen_genes$orthogroup),
                          orthogroup_vec[!orthogroup_vec %in% as.character(Mgen_table$orthogroup)])

Lfla_predictor_genes <- c(as.character(Lfla_genes$orthogroup),
                          orthogroup_vec[!orthogroup_vec %in% as.character(Lfla_table$orthogroup)])



# vector of orthgroups that are predictor in each species, including duplicate                          
common_predictor_genes_vec <- c(Caus_predictor_genes,
                                Ccal_predictor_genes,
                                Pdom_predictor_genes,
                                Pcan_predictor_genes,
                                Mgen_predictor_genes,
                                Lfla_predictor_genes)



# common to all 144 genes now
common_genes <- names(table(common_predictor_genes_vec))[table(common_predictor_genes_vec) == 6]

# plot the ranks for the common predictor genes (lineage colour-coded)
all_species_table %>% filter(orthogroup %in% common_genes) %>% 
  ggplot( aes(x=orthogroup,
              y = nfeatures,
              group = species,
              colour = species)) + 
  geom_point() + scale_colour_manual(values = c('#54278f','#bcbddc', '#fdae6b', 
                                                '#756bb1', '#e6550d','#e6550d')) +
  ggtitle("SVM Ranks of Common Predictor genes") + theme_bw()

# make it into a long format
all_species_table_lf <- all_species_table %>%
  filter(orthogroup %in% common_genes) %>% 
  spread(key = species, value = nfeatures)


# make into matrix (NA are always present in predictor genes so get a rank of 100)
predictor_genes_matrix <- as.matrix(all_species_table_lf[, 2:7])
row.names(predictor_genes_matrix) <- all_species_table_lf$orthogroup
predictor_genes_matrix[is.na(predictor_genes_matrix)] <- 100


# make a heatmap
heatmap(predictor_genes_matrix, scale="none")

# cluster rows
hc.rows <- hclust(dist(predictor_genes_matrix))

plot(hc.rows)

# transpose the matrix and cluster columns
hc.cols <- hclust(dist(t(predictor_genes_matrix)))
plot(hc.cols)
library(RColorBrewer)
coul <- rev(colorRampPalette(brewer.pal(8, "RdPu"))(25))

# draw heatmap with dendrogram of orthologs
# this figure goes in manuscript
pdf("common-predictor-genes-heatmap.pdf")
heatmap(predictor_genes_matrix,
        Colv   = as.dendrogram(hc.cols),
        Rowv   = as.dendrogram(hc.rows),
        scale  = 'none',
        col    = coul, 
        labRow = "")

dev.off()
#############################
### make heatmap with all 3718 orthogroups

# make it into a long format
all_species_3718_table_lf <- all_species_table %>%
  spread(key = species, value = nfeatures)


# make into matrix (NA are always present in predictor genes so get a rank of 100)
all_genes_matrix <- as.matrix(all_species_3718_table_lf[, 2:7])
row.names(all_genes_matrix) <- all_species_3718_table_lf$orthogroup
all_genes_matrix[is.na(all_genes_matrix)] <- 100


# make a heatmap
heatmap(all_genes_matrix, scale="none")

# cluster rows
hc.rows <- hclust(dist(all_genes_matrix))

# colour the signif
plot(hc.rows, colour = "blue")

# transpose the matrix and cluster columns
hc.cols <- hclust(dist(t(all_genes_matrix)))
plot(hc.cols)
#library(RColorBrewer)
coul <- rev(colorRampPalette(brewer.pal(8, "RdPu"))(25))

# make a character vector of colours to mark those predictor genes
colour_predictor_vec <- rep("LightGrey",
                            length(row.names(all_genes_matrix)))
colour_predictor_vec[row.names(all_genes_matrix) %in% 
                       row.names(predictor_genes_matrix)] <- "black"

# draw heatmap with dendrogram of orthologs
# this figure goes in manuscript
pdf("all-predictor-genes-heatmap.pdf")
heatmap(all_genes_matrix,
        Colv   = as.dendrogram(hc.cols),
        Rowv   = as.dendrogram(hc.rows),
        scale  = 'none',
        col    = coul, 
        labRow = "",
        RowSideColors = colour_predictor_vec)

dev.off()


#### BEES only
# draw heatmap with dendrogram of orthologs (only bees)
predictor_genes_matrix_bees <- predictor_genes_matrix[, c(1,2,4)]

heatmap(predictor_genes_matrix_bees,
        Colv = as.dendrogram(hclust(dist(t(predictor_genes_matrix_bees)))),
        Rowv = as.dendrogram(hc.rows),
        scale = 'none')

# plot order of Caus vs order of Ccal
all_species_table %>% 
  filter(species %in% c("Caus", "Ccal") & orthogroup %in% Ccal_genes$orthogroup) %>% 
  spread(key = species, value = nfeatures) %>% 
  ggplot( aes(x=Caus,
              y = Ccal)) + ylim(100, 4000) + xlim(100, 4000) +
  geom_point() + ggtitle("Ranking correlation between the two Ceratinas")

# plot order of Caus vs order of Mgen
all_species_table %>% 
  filter(species %in% c("Caus", "Mgen") & orthogroup %in% Mgen_genes$orthogroup) %>% 
  spread(key = species, value = nfeatures) %>% 
  ggplot( aes(x=Caus,
              y = Mgen)) + ylim(100, 4000) + xlim(100, 4000) +
  geom_point() + ggtitle("Ranking correlation between two bees")


# plot order of Pdom vs order of Pcan
all_species_table %>% 
  filter(species %in% c("Pdom", "Pcan") & orthogroup %in% Pdom_genes$orthogroup) %>% 
  spread(key = species, value = nfeatures) %>% 
  ggplot( aes(x=Pdom,
              y = Pcan)) + ylim(100, 4000) + xlim(100, 4000) +
  geom_point() + ggtitle("Ranking correlation between the two Polistes")

# plot order of Pdom vs order of Lfla
all_species_table %>% 
  filter(species %in% c("Pdom", "Lfla") & orthogroup %in% Pdom_genes$orthogroup) %>% 
  spread(key = species, value = nfeatures) %>% 
  ggplot( aes(x=Pdom,
              y = Lfla)) + ylim(100, 4000) + xlim(100, 4000) +
  geom_point() + ggtitle("Ranking correlation between two wasps")

# plot order of Pdom vs order of Caus
all_species_table %>% 
  filter(species %in% c("Pdom", "Caus") & orthogroup %in% Pdom_genes$orthogroup) %>% 
  spread(key = species, value = nfeatures) %>% 
  ggplot( aes(x=Pdom,
              y = Caus)) + ylim(100, 4000) + xlim(100, 4000) +
  geom_point() + ggtitle("Ranking correlation between lineages")
