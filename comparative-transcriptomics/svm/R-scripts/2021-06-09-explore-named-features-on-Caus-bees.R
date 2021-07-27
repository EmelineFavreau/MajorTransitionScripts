# [Workspace loaded from ~/myriad/projects/MajorTransitionScripts/comparative-transcriptomics/svm/Ceratina_australensisSVM-named_features_bees_only.RData]

# normalised counts from DESeq2
Caus_normalized_counts <- read.delim("~/myriad/projects/MajorTransitionScripts/comparative-transcriptomics/dge-DESeq2/full_set/Ceratina_australensis/Ceratina_australensis_normalized_counts.txt",
                                     stringsAsFactors = FALSE)

# hash table
Caus_3718_gene_orthogroups_list <- read.delim("~/myriad/projects/MajorTransitionScripts/comparative-transcriptomics/svm/input/Ceratina_australensis_3718_gene_orthogroups_list",
                                              header=FALSE, stringsAsFactors = FALSE)

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

# perform analysis
# make a long format dataframe
Caus_long <- tidyr::gather(predictions[complete.cases(predictions),],
                           sample, queen_prob, SRR2917198:SRR2915407,
                           factor_key=FALSE)

# add phenotype column
Caus_long$phenotype <- pheno_data$Phenotype[match(Caus_long$sample,
                                                  pheno_data$SampleName)]

# change any outlier
#Caus_long$queen_prob[Caus_long$queen_prob < 0] <- 0
#Caus_long$queen_prob[Caus_long$queen_prob > 1] <- 1
Caus_long$queen_prob <- as.numeric(Caus_long$queen_prob)
Caus_long$nfeatures <- as.numeric(Caus_long$nfeatures)
# from the output model: number of predictor genes
Caus_predictor_num <- 491

# make a df with means per phenotype
average_pheno <-  Caus_long %>% 
  group_by(phenotype, nfeatures) %>% 
  summarise(queen_prob = mean(queen_prob))

# plot
ggplot(data = Caus_long, 
       aes(x = nfeatures, y = queen_prob, colour = phenotype)) +
  geom_line(aes(group = sample), alpha = .3) + 
  geom_line(data = average_pheno, alpha = .8, size = 3) +
  scale_color_manual(values=c("#6baed6", "#08519c")) +
  scale_x_reverse(limits = c(3718,100)) +
  geom_vline(xintercept = Caus_predictor_num,  col = 'black') +
  ylim(0,1) +
  ggtitle("Ceratina australensis trained on bees only") +
  ylab("Predicted Reproductive Probability") +
  xlab("Number of genes in model") +
  theme_bw()

# the paint drips are not in the same location as the original run
# This means that the outliers are somehow related to the stochasticity of the model


# add gene name to dataframe
Caus_long$gene_name <- as.character(Caus_3718_gene_orthogroups_list$V2[match(Caus_long$orthogroup,
                                                                Caus_3718_gene_orthogroups_list$V1)])
                                    
# add foldchange to dataframe
Caus_long$fold_change <- NA

Caus_long$sample <- as.character(Caus_long$sample)

# iteratively select one gene and one sample, add corresponding fold change
for(this_gene in unique(Caus_long$gene_name)){
    for(this_sample in unique(Caus_long$sample)){
      #print(this_sample)
      #print(this_gene)
      fold_change <- Caus_normalized_counts %>% filter(X == this_gene) %>% 
                      select(this_sample)
      #print(fold_change)
      Caus_long$fold_change[Caus_long$sample == this_sample & 
                              Caus_long$gene_name == this_gene] <- as.vector(unlist(fold_change))
      

    }
}

# 
#Caus_long$sample <- factor(Caus_long$sample)
#Caus_long$fold_change <- order(Caus_long$fold_change)

# Plot input in x, output in y
Caus_long %>%
  ggplot( aes(x = nfeatures,
              y = queen_prob,
              group = sample,
              color = sample,
              shape = phenotype)) +
  geom_point(aes(size = fold_change), alpha = 0.5)


# Plot input in x, output in y
Caus_long %>%
  ggplot( aes(x = fold_change,
              y = queen_prob,
              color = nfeatures,
              shape = sample)) +
  geom_point(alpha = 0.3)

hist(Caus_long$fold_change[Caus_long$queen_prob < 0.25])
hist(Caus_long$fold_change[Caus_long$queen_prob > 0.25])
# their fold change mean are different
t.test(x = Caus_long$fold_change[Caus_long$queen_prob < 0.25],
         y = Caus_long$fold_change[Caus_long$queen_prob > 0.25])
ggplot(data=Caus_long, aes(x=fold_change,
                          group=sample,
                          fill=sample)) +
  geom_density(adjust=1.5, alpha=.4)


unique(Caus_long %>% filter(queen_prob < 0.25) %>% select(orthogroup))
Caus_questionable_genes <- c("OG0004512", "OG0005260", "OG0002904", "OG0003659", "OG0006394")

input_file_for_svm %>% filter(orthogroup %in% Caus_questionable_genes) %>% 
  select(c("SRR2917198_1Aligned.sortedByCoord.out.bam",
           "SRR2917152_1Aligned.sortedByCoord.out.bam",
           "SRR2916699_1Aligned.sortedByCoord.out.bam",
           "SRR2916026_1Aligned.sortedByCoord.out.bam",
           "SRR2916025_1Aligned.sortedByCoord.out.bam",
           "SRR2915407_1Aligned.sortedByCoord.out.bam"))

input_file_for_svm %>%  
  select(c("SRR2917198_1Aligned.sortedByCoord.out.bam",
           "SRR2917152_1Aligned.sortedByCoord.out.bam",
           "SRR2916699_1Aligned.sortedByCoord.out.bam",
           "SRR2916026_1Aligned.sortedByCoord.out.bam",
           "SRR2916025_1Aligned.sortedByCoord.out.bam",
           "SRR2915407_1Aligned.sortedByCoord.out.bam"))

# during leave one out, we performed a variance stabilisation
# check those
svm.counts.vst[Caus_questionable_genes,1:6] %>% colMeans()

svm.counts.vst[,1:6] %>% colMeans()
# I investigated why the paint drips are present
# only 5 orthogroups are responsible, 
# with a sig smaller fold change mean 
# (Welch Two Sample t-test t = -42.098, df = 211.69, p-value < 2.2e-16

# did I add one to an empty orthogroup?
# not

# was it the variance stabilisation at the begining of the svm
# no

# is it a sample size effect?
# maybe. Error rate is positively correlated with ratio test/training sample size ONLY for 1:5 analysis design
# this paper suggest at least 75 samples to start with https://www.sciencedirect.com/science/article/pii/S0003267012016479?via%3Dihub
