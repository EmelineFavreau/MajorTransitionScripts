######################## RANDOMIZATION #####################################


# run SVM on non-randomized data,
# and train on randomized phenotype species 100 times
# keep all error rate for distribution curve

# get libraries
basic_libraries <- c("ggplot2",
                     "tidyverse",
                     "e1071",
                     "DESeq2",
                     "pracma",
                     "limma",
                     "ROCR")

for (lib in basic_libraries) {
  if (require(package = lib, character.only = TRUE)) {
    print("Successful")
  } else {
    print("Installing")
    install.packages(lib)
    library(lib, character.only = TRUE )
  }
}


# load a gene count df with samples in columns and gene id in rownames 
# matrix with orthogroups as row names, samples are column names
# 3718 obs. of  82 variables
counts_clean <- read.csv("../input/readcounts_3718orthogroups_6species.txt",
                         sep = "", stringsAsFactors = FALSE)

# load a df with columns
# Species Phenotype SampleName
pheno_data <- read.delim("../input/species-pheno-sampleName.txt",
                         stringsAsFactors = FALSE)

# replace space by underscore between binomial taxon
pheno_data$Species <- gsub(x = pheno_data$Species,
                           pattern = " ",
                           replacement = "_")

## Aim 1: Define SVM function

## set row names as the orthogroups
row.names(counts_clean) <- counts_clean$orthogroup

# remove orthogroups column
counts_clean <- counts_clean %>% select(-orthogroup)

# sample name should be consistent between datasets
colnames(counts_clean) <- gsub(x           = colnames(counts_clean),
                               pattern     = "_1Aligned.sortedByCoord.out.bam",
                               replacement = "")



# give all 3718 orthogroups
counts_clean_subset <- counts_clean

# change to matrix
counts_clean_subset_mat <- as.matrix(counts_clean_subset)

# set the crossvalidation
crossfold_level  <- 3

# true error rate (Supplementary Table 7)
Caus_error_rate <- 0.1804
Ccal_error_rate <- 0.1917
Mgen_error_rate <- 0.1815
Lfla_error_rate <- 0.1623
Pdom_error_rate <- 0.1798
Pcan_error_rate <- 0.2399

# SET THE SPECIES
focus_species    <- "Liostenogaster_flavolineata"

training_species <- c("Ceratina_calcarata",
  "Polistes_dominula",
  "Megalopta_genalis",
  "Polistes_canadensis",
  "Ceratina_australensis")#,
  #"Liostenogaster_flavolineata")

this_error_rate <- Lfla_error_rate

##############
# make a train data: sample labels and reproductives
# species in training species vec
this_traindata <- pheno_data %>% 
  filter(Species %in% training_species) %>% 
  select(SampleName, Phenotype)

# make a test data: sample labels and reproductives
# just that species
this_testdata <- pheno_data %>% 
  filter(Species == focus_species) %>% 
  select(SampleName, Phenotype)

# instantiate vector for error rate to make a distribution plot
species_error_rate_vec <- rep(0, times = 100)

# test 100 randomized testing datasets
for(i in 1:100){

  # this step is different than the real analysis
  # randomize the samples so they are not the true phenotypes
  this_testdata$SampleName <- sample(x = this_testdata$SampleName,
                                     size = nrow(this_testdata),
                                     replace = FALSE)
  
  # function to train svm
  svm.train <- function(readcounts,
                        traindata,
                        testdata       = NA,
                        referencelevel = "R",
                        kerneltype     = "radial",
                        crossfold      = crossfold_level,
                        vstCheck       = T){
    
    svm.counts.test <- NA
    
    # normalise data
    svm.counts <- readcounts
    
    # perform DESeq's variance stabilizing tranformation, 
    # which is preferable to logging for gene expression data
    if(vstCheck){
      svm.counts.vst <- vst(svm.counts)
    } else {
      svm.counts.vst <- varianceStabilizingTransformation(svm.counts)
    }
    
    # normalize counts between samples 
    svm.counts.vst.quantiles <- normalizeBetweenArrays(svm.counts.vst,
                                                       method = "quantile")
    
    # scale counts and remove zero-variance features
    svm.counts.vst.quantiles.scale <- t(scale(t(svm.counts.vst.quantiles)))
    svm.counts.vst.quantiles.scale <- na.omit(svm.counts.vst.quantiles.scale)
    
    # Divide transcriptomic data into training set 
    # (all other species) and test set (that species) 
    svm.counts.train <- svm.counts.vst.quantiles.scale[ ,
                                                        which(colnames(svm.counts.vst.quantiles.scale) %in% traindata$SampleName)]
    
    if(length(testdata) > 1){
      svm.counts.test <- svm.counts.vst.quantiles.scale[ ,
                                                         which((colnames(svm.counts.vst.quantiles.scale) %in% testdata$SampleName))]
    }
    
    # Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult <- tune("svm", 
                                  train.x     = t(svm.counts.train),
                                  train.y     = as.numeric(traindata$Phenotype == referencelevel),
                                  probability = TRUE, 
                                  scale       = FALSE,
                                  kernel      = kerneltype,
                                  tunecontrol = tune.control(sampling = "cross",
                                                             cross    = crossfold),
                                  ranges = list(gamma = 10^(-7:-5),
                                                cost  = 2^(3:5))
    )
    
    # Final classifier
    svm.counts.classifier <- svm.counts.tuneResult$best.model
    
    svm.counts.prediction <- NULL
    
    if(length(testdata) > 1){
      
      # Make predictions for the test data, if test data were provided.
      svm.counts.prediction <- predict(svm.counts.classifier,
                                       t(svm.counts.test),
                                       type        = "class", 
                                       probability = TRUE)
    }
    
    # output prediction for test data and cross-validation error for training data
    svm.result <- list("prediction"       = svm.counts.prediction,
                       "validation_error" = signif(svm.counts.tuneResult$best.performance, 4),
                       "traincounts"      = svm.counts.train,
                       "testcounts"       = svm.counts.test)
    
    # return results
    return(svm.result)
  }
  
  
  
  
  
  ##  Perform initial classification
  
  
  #### Perform initial classification
  # apply svm to entire set of genes
  svm.full <- svm.train(readcounts  = counts_clean_subset_mat,
                        traindata   = this_traindata, 
                        testdata    = this_testdata, 
                        crossfold   = crossfold_level,
                        vstCheck    = F)
  
  # save error rate to make a distribution plot
  species_error_rate_vec[i] <- svm.full$validation_error
}

# # change to dataset
# data <- data_frame(species_error_rate_vec)
# 
# 
# 
# # plot error rate distribution when randomized testing phenotype
# # with vertical bar for the true error rate
# data %>%
#   ggplot( aes(x=species_error_rate_vec)) +
#   geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
#   geom_vline(xintercept = this_error_rate) +
#   ggtitle(paste("100 SVMs testing randomized", focus_species))
# plotTitle <- paste("../figures/2022-06-16-svm-randomized-",focus_species,".jpg", sep = "")
# ggsave(plotTitle)

# save dataset for quantile calculation
#Liostenogaster_flavolineata_error_rate_vec <- species_error_rate_vec

##########
# for each species, assess the decile in which the true error rate is found
species_error_rate_vec <- Ceratina_australensis_error_rate_vec
# 10%     20%     30%     40%     50%     60%     70%     80%     90% 
# 0.18289 0.18986 0.19517 0.20180 0.20675 0.21328 0.21776 0.22164 0.23208 

species_error_rate_vec <- Ceratina_calcarata_error_rate_vec
# 10%     20%     30%     40%     50%     60%     70%     80%     90% 
# 0.18675 0.19318 0.20125 0.20604 0.21255 0.22082 0.22660 0.23614 0.24521 

species_error_rate_vec <- Megalopta_genalis_error_rate_vec
# 10%     20%     30%     40%     50%     60%     70%     80%     90% 
# 0.18035 0.18850 0.19451 0.20132 0.20915 0.21594 0.22492 0.23282 0.24706 

species_error_rate_vec <- Liostenogaster_flavolineata_error_rate_vec
# 10%     20%     30%     40%     50%     60%     70%     80%     90% 
# 0.18366 0.19556 0.20288 0.21018 0.21790 0.22262 0.23222 0.24086 0.25356

species_error_rate_vec <- Polistes_dominula_error_rate_vec
# 10%     20%     30%     40%     50%     60%     70%     80%     90% 
# 0.18366 0.19556 0.20288 0.21018 0.21790 0.22262 0.23222 0.24086 0.25356

species_error_rate_vec <- Polistes_canadensis_error_rate_vec
# 10%     20%     30%     40%     50%     60%     70%     80%     90% 
# 0.17568 0.18596 0.19278 0.19746 0.20265 0.20914 0.21853 0.22778 0.24925 

# calculate percentile of randomised tests that are worst than the real model
quantile(sort(species_error_rate_vec), probs = seq(.1, .9, by = .1))

# compare with true error rate
# expectation: true error rate should be in the first deciles
# meaning: the true SVM does a better job at predicting than the randomized tests
# Caus_error_rate = 0.1804, in the 1st decile of the error rates distribution from the randomization
# Ccal_error_rate = 0.1917, in the 2nd decile
# Mgen_error_rate = 0.1815, in the 2nd decile
# Lfla_error_rate = 0.1623, in the 1st decile
# Pdom_error_rate = 0.1798, in the 1st decile
# Pcan_error_rate = 0.2399, in the 9th decile
