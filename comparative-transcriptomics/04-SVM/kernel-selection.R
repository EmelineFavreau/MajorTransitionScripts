# checking which kernel to apply to our SVM: linear or radial?
# https://rpubs.com/markloessi/497544

# phenotype predicted based on read count
# focus_species    <- "Ceratina_australensis"
focus_species    <- "Liostenogaster_flavolineata"

  # training_species <- c("Ceratina_calcarata", "Megalopta_genalis")
training_species <- c("Ceratina_calcarata",
                      "Polistes_dominula",
                      "Megalopta_genalis",
                      "Polistes_canadensis",
                      "Ceratina_australensis")#,
                      #"Liostenogaster_flavolineata")


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
pheno_data <- read.delim("input/species-pheno-sampleName.txt",
                         stringsAsFactors = FALSE)

# replace space by underscore between binomial taxon
pheno_data$Species <- gsub(x = pheno_data$Species,
                           pattern = " ",
                           replacement = "_")

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

# scale data
scaled_read_counts <- varianceStabilizingTransformation(counts_clean_subset_mat)

# split the read count in training_set and testing_set
training_set_raw <- scaled_read_counts[, colnames(scaled_read_counts) %in% this_traindata$SampleName]
testing_set_raw <- scaled_read_counts[, colnames(scaled_read_counts) %in% this_testdata$SampleName]

training_set <- as.data.frame(t(training_set_raw))
training_set$Phenotype <- this_traindata$Phenotype[match(row.names(training_set), this_traindata$SampleName)]

testing_set <- as.data.frame(t(testing_set_raw))


# fit svm to training set
classifierR <- svm(formula = Phenotype ~ .,
                  data = training_set,
                  type = 'C-classification', 
                  kernel = 'radial')

# fit svm to training set
classifierL <- svm(formula = Phenotype ~ .,
                   data = training_set,
                   type = 'C-classification', 
                   kernel = 'linear')

# predict the test set result (radial)
y_predR = predict(classifierR, newdata = testing_set)

# predict the test set result (linear)
y_predL = predict(classifierL, newdata = testing_set)

testing_set$Phenotype <- this_testdata$Phenotype[match(row.names(testing_set), this_testdata$SampleName)]
 
# make confusion matrix
cmR <- table(testing_set$Phenotype, y_predR)
cmL <- table(testing_set$Phenotype, y_predL)
cmR
cmL
# calculate accuracy rate
accuracy_rate_radial <- (10+0) / sum(as.numeric(y_predR))
accuracy_rate_radial

# calculate accuracy rate
accuracy_rate_linear <- (10+0) / sum(as.numeric(y_predL))
accuracy_rate_linear

# Pdom: radial (0.5 > 0.4)
# Ceratina australensis: equal
# Ceratina calcarata: equal
# Mgen: radial (0.56 > 0.21)
# Pcan: equal
# Lfl: equal

######################## ROC PLOTS #####################################


# source: https://towardsdatascience.com/roc-curve-and-auc-explained-8ff3438b3154
# source: https://rpubs.com/JanpuHou/359286
# ROC (receiver operating characteristics) curve and 
#AOC (area under the curve) are performance measures 
#that provide a comprehensive evaluation of classification models.



# train the model
x.svm <- svm(formula = Phenotype ~ .,
             data = training_set,
             type = 'C-classification', 
             kernel = 'radial',
             probability = TRUE)

# test the testing data
x.svm.prob <- predict(x.svm,
                      type = "prob",
                      newdata = testing_set,
                      probability = TRUE)

# create prediction objects
x.svm.prob.rocr <- prediction(attr(x.svm.prob,
                                   "probabilities")[,2],
                              this_testdata$Phenotype)

# predictor evaluations are performed: True positive rate AND False positive rate
x.svm.perf <- performance(x.svm.prob.rocr,
                          "tpr","fpr")

# collect info
#Ccal.x.svm.perf <- x.svm.perf
#Pdom.x.svm.perf <- x.svm.perf
#Mgen.x.svm.perf <- x.svm.perf
#Pcan.x.svm.perf <- x.svm.perf
#Caus.x.svm.perf <- x.svm.perf
#Lfl.x.svm.perf <- x.svm.perf

plotName <- "../R-plots/SVM_ROC.pdf"

pdf(plotName)

plot(Caus.x.svm.perf, col=1, main="ROC curves for each species testing")
# Draw a legend.
legend(0,
       1,
       c('C.australensis',
         'C.calcarata',
         'M.genalis',
         'L.flavolineata',
         'P.canadensis',
         'P.dominula'),
       1:6)
plot(Ccal.x.svm.perf, col=2, add=TRUE)
plot(Mgen.x.svm.perf, col=3, add=TRUE)
plot(Lfl.x.svm.perf, col=4, add=TRUE)
plot(Pcan.x.svm.perf, col=5, add=TRUE)
plot(Pdom.x.svm.perf, col=6, add=TRUE)

dev.off()