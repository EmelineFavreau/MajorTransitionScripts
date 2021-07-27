#!/usr/bin/env Rscript
# Support Vector Machine test on one species, by training on other species
# Copyright 2021 Emeline Favreau, University College London.

# Our data are read counts provided by RSEM for 6 species, for 3718 orthogroups.
# this is a categorical SVM (factor), not a regression (numerical)

## Analysis steps:
#- Obtaining data
#- Aim 1: Define SVM function
#- Aim 2: Perform initial classification
#- Aim 3: Perform feature selection
#
args = commandArgs(trailingOnly=TRUE)
# focus_species    <- "Ceratina_australensis"
focus_species    <- args[1]
# training_species <- c("Ceratina_calcarata", "Megalopta_genalis")
training_species <- c(args[-c(1,2)])
# experiment_name <- "bees_only_factors"
experiment_name <- args[2]
crossfold_level  <- 3

# get libraries
basic_libraries <- c("ggplot2",
                     "tidyverse",
                     "e1071",
                     "DESeq2",
                     "pracma",
                     "limma")

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
counts_clean <- read.csv("input/readcounts_3718orthogroups_6species.txt",
                         sep = "", stringsAsFactors = FALSE)

# load a df with columns
# Species Phenotype SampleName
pheno_data <- read.delim("input/species-pheno-sampleName.txt")

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
#counts_clean_subset <- counts_clean[1:1000, ]
# change to matrix
counts_clean_subset_mat <- as.matrix(counts_clean_subset)

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

# create filenames for gene list in full model
GeneListFileName <- paste(focus_species, experiment_name,
                          "predicted_gene_list",
                          sep = "_")
# for plot
plotFileName <- paste(focus_species, experiment_name,
                      "plot.pdf",
                      sep = "_")

# for error rate table
ErrorRateTableFileName <- paste(focus_species, experiment_name,
                                "error_rates",
                                sep = "_")

# for prediction probability table
predictionsTableFileName <- paste(focus_species, experiment_name,
                                  "prediction_table",
                                  sep = "_")

# for workspace
workspaceFileName <- paste(focus_species, experiment_name,
                           ".RData",
                           sep = "")

#### Define SVM function
# the intended inputs for this function are a 2xn dataframe,
# where n is your number of samples. 
# column 1 of the input should be the sample labels in the same order 

# readcounts: 
# an n x m data frame with samples as the columns
# (colnames  = sample names) and genes as the rows

# traindata: 
# a 2 column data frame with column 1 = sample labels for your training data; 
# column 2 = classifications, some of which must match the referencelevel;

# testdata: 
# a 2 column data frame with the same structure as traindata,
# but the first column should contain test samples 
# (the second column is irrelevant because my code is bad)

# kerneltype: svm kernel function, e.g. radial, polynomial, etc.

# crossfold: k for k-fold cross-validation
# which you should vary based on the size of your training sets
# you want to avoid cross-validation bins that contain only one group
# example: training set has 10 samples in 2 classes. 
# If K is 3, we'll have 3 bins of ~3 samples each. 
# If k is 5, 5 bins of 2 samples each. 
# The lower the better if sample size is low.
# I set it at 3 because we're dealing with 3 NR and 3 R

# vstCheck: whether to use vst() or varianceStablizingTransformation();

# usually vst() is fine but it may throw errors for small numbers of genes- 
# then you should set vstCheck to F to use the more stable function


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
                                train.y     = traindata$Phenotype,
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





## Aim 2: Perform initial classification


#### Perform initial classification
# apply svm to entire set of genes
svm.full <- svm.train(readcounts  = counts_clean_subset_mat,
                      traindata   = this_traindata, 
                      testdata    = this_testdata, 
                      crossfold   = crossfold_level,
                      vstCheck    = F)

# check 
print(paste0("Root mean cross-validation error rate for full model: ",
             svm.full$validation_error))



## Aim 3: Perform feature selection

# set the number of error estimates made in each loop in the feature selection 
# (usually 20)
err_esti_num <- 20

# create copy of training data that we can subject to repeated trimming
# while preserving original frame
svm.counts.train.iterate <- svm.full$traincounts
svm.counts.test.iterate <- svm.full$testcounts
# record original number of features
nfeatures <- nrow(svm.counts.train.iterate)

# target number of features 
nfeatures_target <- 100

# set train data
traindata <- this_traindata

# instantiate data frame to hold data on the error of each model
iterations <- data.frame(feature              = character(),
                         error_before_removal = numeric(),
                         stringsAsFactors = FALSE)

# instantiate data frame to hold probability prediction data for each model
predictions <- as.data.frame(matrix(NA,
                                    ncol = length(this_testdata$SampleName)+2))


# name columns
colnames(predictions) <- c("nfeatures", "orthogroup", as.character(this_testdata$SampleName))

# create a list that will contain probabilities for each sample at a given step of feature selection
probabilities_list <- vector(mode = "list", length = nfeatures)

# name each element
names(probabilities_list) <- nfeatures:101

# iteratively remove features until target number is reached
while(nfeatures > nfeatures_target){
  
  error <- c()
  
  #run repeatedly to account for stochasticity in cross-validation
  for(i in 1:err_esti_num){
    #print("starting the loop of 1:1")
    # Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult <- tune("svm", 
                                  train.x     = t(svm.counts.train.iterate), 
                                  train.y     = traindata$Phenotype,
                                  probability = TRUE, 
                                  scale       = FALSE,
                                  kernel      = "radial", 
                                  tunecontrol = tune.control(sampling = "cross", 
                                                             cross  = crossfold_level),
                                  ranges      = list(gamma = 10^(-5:-7),
                                                     cost  = 2^(4:6)))
    
    # record error
    error <- c(error, svm.counts.tuneResult$best.performance)
    
    
    
  }
  
  # sample classifier
  svm.counts.classifier <- svm.counts.tuneResult$best.model
  #print("about to start the predit function")
  
  
  
  # return mean error value
  error <- signif(mean(error), 4)
  
  # extract feature weights
  weights <- (t(svm.counts.classifier$coefs) %*% svm.counts.classifier$SV)
  
  # calculate feature with lowest weight (for ties, choose arbitrarily)
  weakfeature <- colnames(weights)[which(abs(weights) == min(abs(weights)))[1]]
  
  # predict categories of the test data
  this_prediction <- predict(svm.counts.classifier,
          t(svm.counts.test.iterate),
          type = "class", probability = TRUE)
  
  # save all attributes (including probabilities) in the list, including the step and the orthogroup 
  probabilities_list[[nfeatures]] <- c(this_prediction, nfeatures)
  
  # vector of categories for each feature being dropped
  predictions <- rbind(predictions,
                       c(nfeatures, weakfeature,
                         as.vector(this_prediction)))
  
  # remove lowest-weight feature from data frame (train data)
  svm.counts.train.iterate <- subset(svm.counts.train.iterate,
                                     !(rownames(svm.counts.train.iterate) %in%
                                         c(weakfeature)))
  
  # remove lowest-weight feature from data frame (test data)
  svm.counts.test.iterate <- subset(svm.counts.test.iterate,
                                    !(rownames(svm.counts.test.iterate) %in%
                                        c(weakfeature)))
  
  # store removed feature name and error value before removing that feature
  iterations <- rbind(iterations, tibble(feature = weakfeature,
                                         error_before_removal = error))
  # tick down
  nfeatures <- (nfeatures-1)
  #print("finishing the while loop")
  # output every 100 runs to track progress
  if((nfeatures/100)%%1==0){print(paste0("Features remaining: ",
                                         nfeatures))
  }
}

# vector of number of iterations that the model went through
iterLength <- 1:nrow(iterations)

# take moving average to smooth out variation
moving_avg <- movavg(iterations$error_before_removal, 100, "s") 

# plot data to ensure we have the expected 'hockeystick' shape 
hockeyData_plot <- data.frame(num   = iterLength,
                              error = moving_avg)

# x axis is ordered from large to small 
# ie right: when model started; left: when model ended
hockeyData_plot$num <- abs(iterLength - (max(iterLength)+1))

# save pdf of plot
pdf(file = plotFileName)

plot(hockeyData_plot$num, hockeyData_plot$error, xlim = rev(c(0, 1000)))

dev.off()

# get minimum of this curve:
# find the point at which the error window is at its minimum
optimal_removal <- which(moving_avg == min(moving_avg));

# list the features to be removed from the original set of genes
# the set of genes to be removed are ordered from the least important feature
# to the most important feature
features_to_remove <- iterations$feature[1:optimal_removal]

# new dataframe with less-useful features removed
counts_clean_subsample <- subset(counts_clean_subset_mat,
                                 !(rownames(counts_clean_subset_mat) %in%
                                     features_to_remove))

# re-perform support vector classification using the new, 
# optimally caste-separating set of features
svm.optimal <- svm.train(readcounts     = counts_clean_subsample, 
                         referencelevel = "R",
                         traindata      = this_traindata, 
                         testdata       = this_testdata,
                         crossfold      = crossfold_level,
                         vstCheck       = F)

# print important genes in final model
print(paste0("Number of genes included in optimised model: ",
             nrow(counts_clean_subsample)))

# print error rate (to compare with full model error rate printed earlier)
print(paste0("Root mean cross-validation error rate for optimised model: ",
             svm.optimal$validation_error))

# the genes taken in the final model are found rownames(svm.optimal$traincounts)
write.table(x         = svm.optimal$traincounts,
            file      = GeneListFileName,
            row.names = TRUE)     


# save a table with feature to be removed and error rate associated
write.table(x         = iterations,
            file      = ErrorRateTableFileName,
            quote     = FALSE,
            row.names = TRUE)  

# save a table with features to be removed and prob to be a queen for each sample
write.table(x         = predictions,
            file      = predictionsTableFileName,
            quote     = FALSE,
            row.names = TRUE) 

# save workspace; check probabilities_list for more details
save.image(file = workspaceFileName)
