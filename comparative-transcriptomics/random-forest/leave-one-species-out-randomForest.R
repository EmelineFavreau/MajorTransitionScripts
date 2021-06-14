#!/usr/bin/env Rscript
# Random Forest (RF) test on one species, by training on other species
# Copyright 2021 Katherine Geist, Iowa State University; modified from E. Favreau

# Our data are read counts provided by RSEM for 6 species, for 3718 orthogroups.

## Analysis steps:
#- Obtaining data
#- Aim 1: Define RF function
#- Aim 2: Perform initial classification
#- Aim 3: Perform feature selection
#

setwd("/Users/ksg/Dropbox/0.GitHub.Local/MajorTransitionScripts/comparative-transcriptomics/random-forest")

args = commandArgs(trailingOnly=TRUE)
focus_species    <- "Ceratina_australensis"
focus_species    <- args[1]
training_species <- c("Ceratina_calcarata", "Megalopta_genalis")
training_species <- c(args[-c(1,2)])
experiment_name <- "bees_only"
experiment_name <- args[2]
crossfold_level  <- 3

# get libraries
basic_libraries <- c("ggplot2",
                     "tidyverse",
                     "e1071",
                     "randomForest",
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

## Set the seed:
set.seed(seed = 1)   ## Seed would have been 1 in the rf most likely.

# load a gene count df with samples in columns and gene id in rownames 
# matrix with orthogroups as row names, samples are column names
# 3718 obs. of  82 variables
counts_clean <- read.csv("../rf/input/readcounts_3718orthogroups_6species.txt",
                         sep = "", stringsAsFactors = FALSE)

# load a df with columns
# Species Phenotype SampleName
pheno_data <- read.delim("../rf/input/species-pheno-sampleName.txt",
                         stringsAsFactors = FALSE)

# replace space by underscore between binomial taxon
pheno_data$Species <- gsub(x = pheno_data$Species,
                           pattern = " ",
                           replacement = "_")

## Aim 1: Define rf function

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

#### Define Random Forest function
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

?randomForest
# mtry = Number of variables randomly sampled as candidates at each split. 
# Note that the default values are different for classification (sqrt(p) 
# where p is number of variables in x) and regression (p/3)

## Rather than using all predictors and all individuals to make a single tree,
## we will make a forest of many (ntree) trees, each one based on a random
## selection of predictors and individuals.

## Each tree is fit with a bootstrap sample of data (drawn with replacement)
## and "grown" until each node is "pure"

## Each node is split using the best among a subset (of size mtry) of 
## predictors randomly chosen at that node (default is the sqrt of #
## of predictors), but there is a special case called "bagging" that
## uses all predictors.

## The final prediction is made by aggregating across the forest
## by either majority vote or by average.

## Training set = the bootstrap sample, unless we use the bag (full model)
## Test set = "out of bag" or leftover samples
##### out-of-bag error rate is:
########## Permute values of predictor j among all OOB cases
########## Pass OOB data down the tree, save the predictions
########## For case i of OOB and predictor j , get: 
############## OOB error rate with variable j permuted – OOB error rate before permutation
########## Average across forest to get overall variable importance for each predictor j

## Before we mimic the process done for SVM, we want to just know how many our
## genes from the SVM optimized predictions overlap with the bagged RF:
## (So not using the training set at all, just using all the data: counts_clean_subset_mat)
rf <- randomForest(x=t(counts_clean_subset_mat),y=as.factor(pheno_data$Phenotype == "R"),ntree=10000)
## Look at variable importance
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing=TRUE)
plot(c(1:nrow(counts_clean_subset_mat)),imp.temp[t],log='x',cex.main=1.5,
     xlab='gene rank',ylab='variable importance',cex.lab=1.5,
     pch=16,main='ALL subset results')



rf.train <- function(readcounts,
                      traindata,
                      testdata       = NA,
                      referencelevel = "R",
                      crossfold      = crossfold_level,
                      vstCheck       = T){
  
  rf.counts.test <- NA
  
  # normalise data
  rf.counts <- readcounts
  
  # perform DESeq's variance stabilizing tranformation, 
  # which is preferable to logging for gene expression data
  if(vstCheck){
    rf.counts.vst <- vst(rf.counts)
  } else {
    rf.counts.vst <- varianceStabilizingTransformation(rf.counts)
  }
  
  # normalize counts between samples 
  rf.counts.vst.quantiles <- normalizeBetweenArrays(rf.counts.vst,
                                                     method = "quantile")
  
  # scale counts and remove zero-variance features
  rf.counts.vst.quantiles.scale <- t(scale(t(rf.counts.vst.quantiles)))
  rf.counts.vst.quantiles.scale <- na.omit(rf.counts.vst.quantiles.scale)
  
  # Divide transcriptomic data into training set 
  # (all other species) and test set (that species) 
  rf.counts.train <- rf.counts.vst.quantiles.scale[ ,
      which(colnames(rf.counts.vst.quantiles.scale) %in% traindata$SampleName)]
  
  if(length(testdata) > 1){
    rf.counts.test <- rf.counts.vst.quantiles.scale[ ,
     which((colnames(rf.counts.vst.quantiles.scale) %in% testdata$SampleName))]
  }
  
  # Perform a grid search to optimise SVM parameters
  rf.counts.tuneResult <- tune("svm", 
                train.x     = t(rf.counts.train),
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
  rf.counts.classifier <- rf.counts.tuneResult$best.model
  
  rf.counts.prediction <- NULL
  
  if(length(testdata) > 1){
    
    # Make predictions for the test data, if test data were provided.
    rf.counts.prediction <- predict(rf.counts.classifier,
                                     t(rf.counts.test),
                                     type        = "class", 
                                     probability = TRUE)
  }
  
  # output prediction for test data and cross-validation error for training data
  rf.result <- list("prediction"       = rf.counts.prediction,
                     "validation_error" = signif(rf.counts.tuneResult$best.performance, 4),
                     "traincounts"      = rf.counts.train,
                     "testcounts"       = rf.counts.test)
  
  # return results
  return(rf.result)
}





## Aim 2: Perform initial classification


#### Perform initial classification
# apply rf to entire set of genes
rf.full <- rf.train(readcounts  = counts_clean_subset_mat,
                      traindata   = this_traindata, 
                      testdata    = this_testdata, 
                      crossfold   = crossfold_level,
                      vstCheck    = F)

# check 
print(paste0("Root mean cross-validation error rate for full model: ",
             rf.full$validation_error))



## Aim 3: Perform feature selection

# set the number of error estimates made in each loop in the feature selection 
# (usually 20)
err_esti_num <- 20

# create copy of training data that we can subject to repeated trimming
# while preserving original frame
rf.counts.train.iterate <- rf.full$traincounts
rf.counts.test.iterate <- rf.full$testcounts
# record original number of features
nfeatures <- nrow(rf.counts.train.iterate)

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
                                   ncol = length(this_testdata$SampleName)+1))
 
# name columns
colnames(predictions) <- c("nfeatures", this_testdata$SampleName)


# iteratively remove features until target number is reached
while(nfeatures > nfeatures_target){
  
  error <- c()
  
  #run repeatedly to account for stochasticity in cross-validation
  for(i in 1:err_esti_num){
    #print("starting the loop of 1:1")
    # Perform a grid search to optimise rf parameters
    rf.counts.tuneResult <- tune("rf", 
                          train.x     = t(rf.counts.train.iterate), 
                          train.y     = as.numeric(traindata$Phenotype == "R"),
                          probability = TRUE, 
                          scale       = FALSE,
                          kernel      = "radial", 
                          tunecontrol = tune.control(sampling = "cross", 
                                                     cross  = crossfold_level),
                          ranges      = list(gamma = 10^(-5:-7),
                                             cost  = 2^(4:6)))
    
    # record error
    error <- c(error, rf.counts.tuneResult$best.performance)
    
    
    
  }
  
  # sample classifier
  rf.counts.classifier <- rf.counts.tuneResult$best.model
  #print("about to start the predit function")
  # vector of probability of being reproductive
  predictions <- rbind(predictions,
                       c(nfeatures,
                         predict(rf.counts.classifier,
                  t(rf.counts.test.iterate), type = "class", probability = TRUE)))
  
  # return mean error value
  error <- signif(mean(error), 4)
  
  # extract feature weights
  weights <- (t(rf.counts.classifier$coefs) %*% rf.counts.classifier$SV)
  
  # calculate feature with lowest weight (for ties, choose arbitrarily)
  weakfeature <- colnames(weights)[which(abs(weights) == min(abs(weights)))[1]]
  
  # remove lowest-weight feature from data frame (train data)
  rf.counts.train.iterate <- subset(rf.counts.train.iterate,
                                     !(rownames(rf.counts.train.iterate) %in%
                                         c(weakfeature)))
  
  # remove lowest-weight feature from data frame (test data)
  rf.counts.test.iterate <- subset(rf.counts.test.iterate,
                                     !(rownames(rf.counts.test.iterate) %in%
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
rf.optimal <- rf.train(readcounts     = counts_clean_subsample, 
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
             rf.optimal$validation_error))

# the genes taken in the final model are found rownames(rf.optimal$traincounts)
write.table(x         = rf.optimal$traincounts,
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

# save workspace
save.image(file = workspaceFileName)
