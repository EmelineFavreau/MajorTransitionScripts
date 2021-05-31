#---
#title: "SVM test on Megalopta genalis, by training on bee only"
#author: "Emeline Favreau"
#date: "2021/05/15"

#---
#### Copyright 2021 Emeline Favreau, University College London.

#Based on https://github.com/BenjaminATaylor/Taylor-et-al-2020-demo
#https://cran.r-project.org/web/packages/e1071/vignettes/svmdoc.pdf
#https://towardsdatascience.com/a-guide-to-svm-parameter-tuning-8bfe6b8a452c
#https://scikit-learn.org/stable/modules/cross_validation.html
#https://rstudio-pubs-static.s3.amazonaws.com/271792_96b51b7fa2af4b3f808d04f3f3051516.html
#https://towardsdatascience.com/cross-validation-430d9a5fee22



#---
## Objective of analysis
#This is a test script to test the code provided by Ben.
#Our data are read counts provided by RSEM for 6 species, for 3718 orthogroups.
#We use one species as test data, and the three bee species as training data.


## Analysis steps:
#- Obtaining data
#- Aim 1: Define SVM function
#- Aim 2: Perform initial classification
#- Aim 3: Perform feature selection




#```{r load all the libraries, eval = TRUE, echo = FALSE, include = FALSE}
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

#BiocManager::install("limma")
#```



#```{r import data, eval = TRUE, echo = FALSE, include = FALSE}
# load a gene count df with samples in columns and gene id in rownames (1000 orthogroups, 96 samples)

# matrix with orthogroups as row names, samples are column names
# 3718 obs. of  82 variables
counts_clean <- read.csv("input/readcounts_3718orthogroups_6species.txt",
                         sep ="", stringsAsFactors = FALSE)

# load a df with columns
# Species Phenotype SampleName
pheno_data <- read.delim("input/species-pheno-sampleName.txt", stringsAsFactors = FALSE)
#```

## Aim 1: Define SVM function

#```{r aim 1 define svm function, eval = TRUE, echo = FALSE, include = TRUE}

## as a toydataset, subset the data to 1000 orthogroups to start
row.names(counts_clean) <- counts_clean$orthogroup

counts_clean <- counts_clean %>% select(-orthogroup)

# sample name should be consistent between datasets
colnames(counts_clean) <- gsub(x = colnames(counts_clean),
                               pattern = "_1Aligned.sortedByCoord.out.bam",
                               replacement = "")

# subsample to 1000 orthogroups
#counts_clean_subset <- counts_clean[1:1000, ]

#alternatively give all 3718 orthogroups
counts_clean_subset <- counts_clean

# change to matrix
counts_clean_subset_mat <- as.matrix(counts_clean_subset)

# check order
# colnames(counts_clean) %in% pheno_data$SampleName

# make a train data: sample labels and reproductives
# all bee species
this_traindata <- pheno_data %>% 
  filter(Species %in% c("Ceratina australensis", "Ceratina calcarata")) %>% 
  select(SampleName, Phenotype)

# make a test data: sample labels and reproductives
# just calcarata
this_testdata <- pheno_data %>% 
  filter(Species == "Megalopta genalis") %>% 
  select(SampleName, Phenotype)

# create filename for gene list in full model
FileName <- "Megalopta_genalis_svm_agst_bee_only_predicted_gene_list"
plotFileName <- "Megalopta_genalis_svm_agst_bee_only_plot.pdf"
ErrorRateTable <- "Megalopta_genalis_svm_agst_bee_only_error_rates"
#### Define SVM function
# the intended inputs for this function are a 2xn dataframe,
# where n is your number of samples. 
# column 1 of the input should be the sample labels *in the same order as they appear as 

# readcounts: an n x m data frame with samples as the columns (colnames  = sample names) and genes as the rows

# traindata: a 2 column data frame with column 1 = sample labels for your training data; 
# column 2 = classifications, some of which must match the referencelevel;

# testdata: a 2 column data frame with the same structure as traindata, but the first column should contain test samples (the second column is irrelevant because my code is bad)

# kerneltype: svm kernel function, e.g. radial, polynomial, etc.; I find radial is reliable 

# crossfold: k for k-fold cross-validation
# which you should vary based on the size of your training sets
# you want to avoid cross-validation bins that contain only one group
#example: training set has 10 samples in 2 classes. If K is 3, we'll have 3 bins of ~3 samples each. If k is 5, 5 bins of 2 samples each. The lower the better if sample size is low.
# I set it at 2 because we're dealing with 3 NR and 3 R

# vstCheck: whether to use vst() or varianceStablizingTransformation();

# usually vst() is fine but it may throw errors for small numbers of genes- 
# then you should set vstCheck to F to use the more stable function


svm.train <- function(readcounts,
                      traindata,
                      testdata       = NA,
                      referencelevel = "R",
                      kerneltype     = "radial",
                      crossfold      = 2,
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
  # (queens and workers from control) and test set (individuals from treatment) 
  svm.counts.train <- svm.counts.vst.quantiles.scale[, which(colnames(svm.counts.vst.quantiles.scale) %in% traindata$SampleName)]
  
  if(length(testdata) > 1){
    svm.counts.test <- svm.counts.vst.quantiles.scale[ ,which((colnames(svm.counts.vst.quantiles.scale) %in% testdata$SampleName))]
  }
  
  # Perform a grid search to optimise SVM parameters
  svm.counts.tuneResult <- tune("svm", 
                               train.x     = t(svm.counts.train),
                               train.y     = as.numeric(traindata$Phenotype == referencelevel),
                               probability = TRUE, 
                               scale       = FALSE,
                               kernel      = kerneltype,
                               tunecontrol = tune.control(sampling = "cross",
                                                          cross = crossfold),
                               ranges = list(gamma = 10^(-7:-5),
                                             cost = 2^(3:5))
  )
  
  # Final classifier
  svm.counts.classifier <- svm.counts.tuneResult$best.model
  
  svm.counts.prediction <- NULL
  
  if(length(testdata) > 1){
    
    # Make predictions for the test data, if test data were provided.
    svm.counts.prediction <- predict(svm.counts.classifier,
                                    t(svm.counts.test),
                                    type = "class", 
                                    probability = TRUE)
  }
  
  # output prediction for test data and cross-validation error for training data
  svm.result <- list("prediction"      = svm.counts.prediction,
                    "validation_error" = signif(svm.counts.tuneResult$best.performance, 4),
                    "traincounts"      = svm.counts.train,
                    "testcounts"       = svm.counts.test)
  
  # return results
  return(svm.result)
}

#```



## Aim 2: Perform initial classification

#```{r aim 2 Perform initial classification, eval = TRUE, echo = FALSE, include = TRUE}
#### Perform initial classification
# apply svm to entire set of genes
svm.full <- svm.train(readcounts  = counts_clean_subset_mat,
                      traindata   = this_traindata, 
                      testdata    = this_testdata, 
                      crossfold   = 2,
                      vstCheck    = F)

# check 
print(paste0("Root mean cross-validation error rate for full model: ",
             svm.full$validation_error))
#```
# Root mean cross-validation error rate for full model: 0.2635
 
# We aim for the lowest root mean cross-validation error rate for full model.



## Aim 3: Perform feature selection

#```{r aim 3 Perform feature selection, eval = TRUE, echo = FALSE, include = TRUE}
#### Perform feature selection

# for the purpose of testing, 
# set the number of error estimates made in each loop in the feature selection 
# (usually 20)
err_esti_num <- 20

# create copy of training data that we can subject to repeated trimming
# while preserving original frame
svm.counts.train.iterate <- svm.full$traincounts

# record original number of features
nfeatures <- nrow(svm.counts.train.iterate)

# target number of features 
nfeatures_target <- 100

# set train data
traindata <- this_traindata

# instantiate data frame to hold data on the error of each model
iterations <- data.frame(feature = character(),
                        error_before_removal = numeric())

# iteratively remove features until target number is reached
while(nfeatures > nfeatures_target){
  
  error <- c()
  
  #run repeatedly to account for stochasticity in cross-validation
  for(i in 1:err_esti_num){
    
    # Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult <- tune("svm", 
                                 train.x     = t(svm.counts.train.iterate), 
                                 train.y     =  as.numeric(traindata$Phenotype == "R"),
                                 probability = TRUE, 
                                 scale       = FALSE,
                                 kernel      = "radial", 
                                 tunecontrol = tune.control(sampling = "cross", 
                                                            cross = 2),
                                 ranges = list(gamma = 10^(-5:-7),
                                               cost = 2^(4:6)))
    
    # record error
    error <- c(error, svm.counts.tuneResult$best.performance)
  }
  
  # sample classifier
  svm.counts.classifier <- svm.counts.tuneResult$best.model
  
  # return mean error value
  error <- signif(mean(error), 4)
  
  # extract feature weights
  weights <- (t(svm.counts.classifier$coefs) %*% svm.counts.classifier$SV)
  
  # calculate feature with lowest weight (for ties, choose arbitrarily)
  weakfeature <- colnames(weights)[which(abs(weights) == min(abs(weights)))[1]]
  
  # remove lowest-weight feature from data frame
  svm.counts.train.iterate <- subset(svm.counts.train.iterate,
                                     !(rownames(svm.counts.train.iterate) %in%
                                         c(weakfeature)))
  
  # in a dataframe, store removed feature name and error value before removing that feature
  iterations <- rbind(iterations, tibble(feature = weakfeature,
                                        error_before_removal = error))
  # tick down
  nfeatures <- (nfeatures-1)
  
  # output every 100 runs to track progress
  if((nfeatures/100)%%1==0){print(paste0("Features remaining: ",
                                         nfeatures))}
}

iterLength <- 1:nrow(iterations)

# take moving average to smooth out variation
moving_avg <- movavg(iterations$error_before_removal, 100, "s") 

# plot data to ensure we have the expected 'hockeystick' shape 
# note that this will be truncated for the subsetted dataset provided for demoing!
hockeyData <- data.frame(num = iterLength,
                        error = moving_avg)
hockeyData_plot <- hockeyData

hockeyData_plot$num <- abs(iterLength - (max(iterLength)+1))

# save pdf of plot
pdf(file=plotFileName)

plot(hockeyData_plot$num, hockeyData_plot$error, xlim = rev(c(0, 1000)))

dev.off()

# get minimum of this curve to find the point at which the error window is at its minimum
optimal_removal <- which(moving_avg == min(moving_avg));

# list the features to be removed from the original set of genes
features_to_remove <- iterations$feature[1:optimal_removal]

# new dataframe with less-useful features removed
counts_clean_subsample <- subset(counts_clean_subset_mat,
                !(rownames(counts_clean_subset_mat) %in% features_to_remove))

# re-perform support vector classification using the new, 
# optimally caste-separating set of features
svm.optimal <- svm.train(counts_clean_subsample, 
                        referencelevel = "R",
                        this_traindata, 
                        this_testdata,
                        crossfold = 2,
                        vstCheck = F)

# 
print(paste0("Number of genes included in optimised model: ",
             nrow(counts_clean_subsample)))

print(paste0("Root mean cross-validation error rate for optimised model: ",
             svm.optimal$validation_error))

# the genes taken in the final model are found rownames(svm.optimal$traincounts)
write.table(x = svm.optimal$traincounts,
            file = FileName,
            row.names = TRUE)     

# save a table with feature to be removed and error rate associated
write.table(x = iterations,
            file = ErrorRateTable,
            quotes = FALSE,
            row.names = TRUE)  

