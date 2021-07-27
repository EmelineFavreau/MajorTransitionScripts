# script to filter read count files ready for SVM analysis
# input: one file per species, with all samples
# output: one file, subsetted to orthogroups and samples needed
# other useful input files: list of orthogrop/specie sgene hash tables, sample and phenotype info


library("tidyverse")

# load all species' read count data
Ceratina_australensis_rawReadCounts <- read.table("input/Ceratina-australensis_raw-read-counts",
                       header = TRUE, stringsAsFactors = FALSE)

Ceratina_calcarata_rawReadCounts <- read.table("input/Ceratina-calcarata_raw-read-counts",
                       header = TRUE, stringsAsFactors = FALSE)

Megalopta_genalis_rawReadCounts <- read.table("input/Megalopta-genalis_raw-read-counts",
                                              header = TRUE, stringsAsFactors = FALSE)
Liostenogaster_flavolineata_rawReadCounts <- read.table("input/Liostenogaster-flavolineatas_raw-read-counts",
header = TRUE, stringsAsFactors = FALSE)

Polistes_canadensis_rawReadCounts <- read.table("input/Polistes-canadensis_raw-read-counts",
header = TRUE, stringsAsFactors = FALSE)

Polistes_dominula_rawReadCounts <- read.table("input/Polistes-dominula_raw-read-counts",
                       header = TRUE, stringsAsFactors = FALSE)

# load all species' orthogroup data
Ceratina_australensis_orthogroupsList <- read.table("input/Ceratina_australensis_3718_gene_orthogroups_list",
              header = FALSE, stringsAsFactors = FALSE)

Ceratina_calcarata_orthogroupsList <- read.table("input/Ceratina_calcarata_3718_gene_orthogroups_list",
              header = FALSE, stringsAsFactors = FALSE)

Megalopta_genalis_orthogroupsList <- read.table("input/Megalopta_genalis_3718_gene_orthogroups_list",
                                                header = FALSE, stringsAsFactors = FALSE)
Liostenogaster_flavolineata_orthogroupsList <- read.table("input/Liostenogaster_flavolineata_3718_gene_orthogroups_list",
header = FALSE, stringsAsFactors = FALSE)

Polistes_canadensis_orthogroupsList <- read.table("input/Polistes_canadensis_3718_gene_orthogroups_list",
header = FALSE, stringsAsFactors = FALSE)

Polistes_dominula_orthogroupsList <- read.table("input/Polistes_dominula_3718_gene_orthogroups_list",
              header = FALSE, stringsAsFactors = FALSE)

# load species names, reproductive status and sample name
speciesPhenoSampleName <- read.delim("input/species-pheno-sampleName.txt",
                                     stringsAsFactors = FALSE)


# update sample names for this species (this will be the columns of input data to svm)
speciesPhenoSampleName$SampleName <- gsub(pattern = "$",
                                          replacement = "_1Aligned.sortedByCoord.out.bam",
                                          x = speciesPhenoSampleName$SampleName)

# make a input data to svm
input_file_for_svm <- as.data.frame(matrix(0,
                                           nrow = length(Ceratina_australensis_orthogroupsList$V1),
                                           ncol = (1+nrow(speciesPhenoSampleName))),
                                    stringsAsFactors = FALSE)

# name the columns
colnames(input_file_for_svm) <- c("orthogroup",
                                  speciesPhenoSampleName$SampleName)

# add the orthogroups
input_file_for_svm$orthogroup <- Ceratina_australensis_orthogroupsList$V1
  
# vector of species
species_vec <- unique(speciesPhenoSampleName$Species)

###############################################################################
# fill the input with read counts based on species
# australensis
species          <- "Ceratina australensis"
rawReadCounts    <- Ceratina_australensis_rawReadCounts
orthogroupsList  <- Ceratina_australensis_orthogroupsList
  
  
# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName %>%
    filter(Species == species) 
  
# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                  all_of(speciesSamplesdf$SampleName)) %>% 
    filter(Geneid %in% orthogroupsList$V2)
  
# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                             orthogroupsList$V2)]
  
# add australensis info into the svm input file
  for(this_orthogroup in subsetrawReadCounts$orthogroup){
    for(this_sample in speciesSamplesdf$SampleName){
      input_file_for_svm[input_file_for_svm$orthogroup == this_orthogroup,
                              colnames(input_file_for_svm) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                                  colnames(subsetrawReadCounts) == this_sample]
    }
  } 



#####################
# calcarata
species          <- "Ceratina calcarata"
rawReadCounts    <- Ceratina_calcarata_rawReadCounts
orthogroupsList  <- Ceratina_calcarata_orthogroupsList


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm[input_file_for_svm$orthogroup == this_orthogroup,
                       colnames(input_file_for_svm) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                           colnames(subsetrawReadCounts) == this_sample]
  }
} 


#####################
# genalis
species          <- "Megalopta genalis"
rawReadCounts    <- Megalopta_genalis_rawReadCounts
orthogroupsList  <- Megalopta_genalis_orthogroupsList


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm[input_file_for_svm$orthogroup == this_orthogroup,
                       colnames(input_file_for_svm) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                           colnames(subsetrawReadCounts) == this_sample]
  }
} 


#####################
# Liostenogaster flavolineata
species          <- "Liostenogaster flavolineata"
rawReadCounts    <- Liostenogaster_flavolineata_rawReadCounts
orthogroupsList  <- Liostenogaster_flavolineata_orthogroupsList


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm[input_file_for_svm$orthogroup == this_orthogroup,
                       colnames(input_file_for_svm) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                           colnames(subsetrawReadCounts) == this_sample]
  }
} 

#####################
# Polistes canadensis
species          <- "Polistes canadensis"
rawReadCounts    <- Polistes_canadensis_rawReadCounts
orthogroupsList  <- Polistes_canadensis_orthogroupsList


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm[input_file_for_svm$orthogroup == this_orthogroup,
                       colnames(input_file_for_svm) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                           colnames(subsetrawReadCounts) == this_sample]
  }
} 

#####################
# Polistes dominula
species          <- "Polistes dominula"
rawReadCounts    <- Polistes_dominula_rawReadCounts
orthogroupsList  <- Polistes_dominula_orthogroupsList


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm[input_file_for_svm$orthogroup == this_orthogroup,
                       colnames(input_file_for_svm) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                           colnames(subsetrawReadCounts) == this_sample]
  }
} 

############
# combine the first colum (orthogroup) with the following result
# add 1 to all values, so that there are no zero counts
input_file_for_svm1 <- cbind(input_file_for_svm[, 1],
                            input_file_for_svm[, -1] + 1)
# name columns
colnames(input_file_for_svm1) <- colnames(input_file_for_svm)

# check result                
# str(input_file_for_svm)
write.table(file = "input/readcounts_3718orthogroups_6species.txt",
            x = input_file_for_svm, row.names = FALSE)

################
################
################
# Make the input file for the wasp-only orthogroups

# load all species' orthogroup data

Liostenogaster_flavolineata_orthogroupsList6983 <- read.table("input/Liostenogaster_flavolineata_6983_gene_orthogroups_list",
                                                              header = FALSE, stringsAsFactors = FALSE)

Polistes_canadensis_orthogroupsList6983 <- read.table("input/Polistes_canadensis_6983_gene_orthogroups_list",
                                                      header = FALSE, stringsAsFactors = FALSE)

Polistes_dominula_orthogroupsList6983 <- read.table("input/Polistes_dominula_6983_gene_orthogroups_list",
                                                    header = FALSE, stringsAsFactors = FALSE)

# update file for just the wasps
speciesPhenoSampleName6983 <- speciesPhenoSampleName %>% filter(Species %in% c("Polistes dominula",
                                                                               "Polistes canadensis",
                                                                               "Liostenogaster flavolineata"))

# make a input data to svm
input_file_for_svm6983 <- as.data.frame(matrix(0,
                                               nrow = length(Polistes_canadensis_orthogroupsList6983$V1),
                                               ncol = (1+nrow(speciesPhenoSampleName6983))),
                                        stringsAsFactors = FALSE)

# name the columns
colnames(input_file_for_svm6983) <- c("orthogroup",
                                      speciesPhenoSampleName6983$SampleName)

# add the orthogroups
input_file_for_svm6983$orthogroup <- Polistes_canadensis_orthogroupsList6983$V1

# vector of species
species_vec <- unique(speciesPhenoSampleName6983$Species)

###############################################################################
# fill the input with read counts based on species


#####################
# Liostenogaster flavolineata
species          <- "Liostenogaster flavolineata"
rawReadCounts    <- Liostenogaster_flavolineata_rawReadCounts
orthogroupsList  <- Liostenogaster_flavolineata_orthogroupsList6983


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName6983 %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm6983[input_file_for_svm6983$orthogroup == this_orthogroup,
                           colnames(input_file_for_svm6983) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                                   colnames(subsetrawReadCounts) == this_sample]
  }
} 

#####################
# Polistes canadensis
species          <- "Polistes canadensis"
rawReadCounts    <- Polistes_canadensis_rawReadCounts
orthogroupsList  <- Polistes_canadensis_orthogroupsList6983


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName6983 %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm6983[input_file_for_svm6983$orthogroup == this_orthogroup,
                           colnames(input_file_for_svm6983) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                                   colnames(subsetrawReadCounts) == this_sample]
  }
} 

#####################
# Polistes dominula
species          <- "Polistes dominula"
rawReadCounts    <- Polistes_dominula_rawReadCounts
orthogroupsList  <- Polistes_dominula_orthogroupsList6983


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName6983 %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm6983[input_file_for_svm6983$orthogroup == this_orthogroup,
                           colnames(input_file_for_svm6983) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                                   colnames(subsetrawReadCounts) == this_sample]
  }
} 

############
# combine the first colum (orthogroup) with the following result
# add 1 to all values, so that there are no zero counts
input_file_for_svm69831 <- cbind(input_file_for_svm6983[, 1],
                                 input_file_for_svm6983[, -1] + 1)
# name columns
colnames(input_file_for_svm69831) <- colnames(input_file_for_svm6983)

# check result                
# str(input_file_for_svm69831)
write.table(file = "input/readcounts_6983orthogroups_wasponly.txt",
            x = input_file_for_svm69831, row.names = FALSE)



################
################
################
# Make the input file for the bee-only orthogroups

# load all species' orthogroup data

Megalopta_genalis_orthogroupsList5787 <- read.table("input/Megalopta_genalis_5787_gene_orthogroups_list",
                                                              header = FALSE, stringsAsFactors = FALSE)

Ceratina_australensis_orthogroupsList5787 <- read.table("input/Ceratina_australensis_5787_gene_orthogroups_list",
                                                      header = FALSE, stringsAsFactors = FALSE)

Ceratina_calcarata_orthogroupsList5787 <- read.table("input/Ceratina_calcarata_5787_gene_orthogroups_list",
                                                    header = FALSE, stringsAsFactors = FALSE)

# update file for just the bees
speciesPhenoSampleName5787 <- speciesPhenoSampleName %>% filter(Species %in% c("Ceratina calcarata",
                                                                               "Ceratina australensis",
                                                                               "Megalopta genalis"))

# make a input data to svm
input_file_for_svm5787 <- as.data.frame(matrix(0,
                                               nrow = length(Ceratina_australensis_orthogroupsList5787$V1),
                                               ncol = (1+nrow(speciesPhenoSampleName5787))),
                                        stringsAsFactors = FALSE)

# name the columns
colnames(input_file_for_svm5787) <- c("orthogroup",
                                      speciesPhenoSampleName5787$SampleName)

# add the orthogroups
input_file_for_svm5787$orthogroup <- Ceratina_australensis_orthogroupsList5787$V1

# vector of species
species_vec <- unique(speciesPhenoSampleName5787$Species)

###############################################################################
# fill the input with read counts based on species


#####################
# Megalopta genalis
species          <- "Megalopta genalis"
rawReadCounts    <- Megalopta_genalis_rawReadCounts
orthogroupsList  <- Megalopta_genalis_orthogroupsList5787


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName5787 %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm5787[input_file_for_svm5787$orthogroup == this_orthogroup,
                           colnames(input_file_for_svm5787) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                                   colnames(subsetrawReadCounts) == this_sample]
  }
} 

#####################
# Ceratina australensis
species          <- "Ceratina australensis"
rawReadCounts    <- Ceratina_australensis_rawReadCounts
orthogroupsList  <- Ceratina_australensis_orthogroupsList5787


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName5787 %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm5787[input_file_for_svm5787$orthogroup == this_orthogroup,
                           colnames(input_file_for_svm5787) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                                   colnames(subsetrawReadCounts) == this_sample]
  }
} 

#####################
# Ceratina calcarata
species          <- "Ceratina calcarata"
rawReadCounts    <- Ceratina_calcarata_rawReadCounts
orthogroupsList  <- Ceratina_calcarata_orthogroupsList5787


# obtain sample info for this species
speciesSamplesdf <- speciesPhenoSampleName5787 %>%
  filter(Species == species) 

# subset raw read counts for these samples and for only the orthogroups
subsetrawReadCounts <- rawReadCounts %>% select(Geneid,
                                                speciesSamplesdf$SampleName) %>% 
  filter(Geneid %in% orthogroupsList$V2)

# add the orthogroups column
subsetrawReadCounts$orthogroup <- orthogroupsList$V1[match(subsetrawReadCounts$Geneid, 
                                                           orthogroupsList$V2)]

# add species info into the svm input file
for(this_orthogroup in subsetrawReadCounts$orthogroup){
  for(this_sample in speciesSamplesdf$SampleName){
    input_file_for_svm5787[input_file_for_svm5787$orthogroup == this_orthogroup,
                           colnames(input_file_for_svm5787) == this_sample] <- subsetrawReadCounts[subsetrawReadCounts$orthogroup == this_orthogroup,
                                                                                                   colnames(subsetrawReadCounts) == this_sample]
  }
} 

############
# combine the first colum (orthogroup) with the following result
# add 1 to all values, so that there are no zero counts
input_file_for_svm57871 <- cbind(input_file_for_svm5787[, 1],
                                 input_file_for_svm5787[, -1] + 1)
# name columns
colnames(input_file_for_svm57871) <- colnames(input_file_for_svm5787)

# check result                
# str(input_file_for_svm57871)
write.table(file = "input/readcounts_5787orthogroups_beeonly.txt",
            x = input_file_for_svm57871, row.names = FALSE)
