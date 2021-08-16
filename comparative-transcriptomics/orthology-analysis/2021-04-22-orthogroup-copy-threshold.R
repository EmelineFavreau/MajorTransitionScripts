######################
#install.packages("UpSetR")
library(UpSetR)
#https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
library("tidyverse")
library(ggplot2)

# import the data for all Species (including drosophila and apis)
Orthogroups.GeneCount <- read.delim("~/Orthogroups.GeneCount.tsv")

# select rows if values are between 1 and 3
filtered_mulitple_copies_df <- Orthogroups.GeneCount %>% filter(Apis_mellifera.longest.isoforms 
                                 %in% c(1, 2, 3))  %>%     
  filter(Ceratina_australensis.longest.isoforms %in% c(1, 2, 3))  %>%        
  filter(Ceratina_calcarata.longest.isoforms %in% c(1, 2, 3))  %>%            
  filter(Drosophila_melanogaster %in% c(1, 2, 3))  %>%                      
  filter(Liostenogaster_flavolineata.longest.isoforms %in% c(1, 2, 3))  %>%   
  filter(Megalopta_genalis.longest.isoforms %in% c(1, 2, 3))  %>%          
  filter(Polistes_canadensis.longest.isoforms %in% c(1, 2, 3))  %>%         
  filter(Polistes_dominula.longest.isoforms %in% c(1, 2, 3))
nrow(filtered_mulitple_copies_df) # 948

# select rows if values are between 0 or 1
filtered_maybe_singlecopy_df <- Orthogroups.GeneCount %>% filter(Apis_mellifera.longest.isoforms 
                                                %in% c(0, 1))  %>%     
  filter(Ceratina_australensis.longest.isoforms %in% c(0, 1))  %>%        
  filter(Ceratina_calcarata.longest.isoforms %in% c(0, 1))  %>%            
  filter(Drosophila_melanogaster %in% c(0, 1))  %>%                      
  filter(Liostenogaster_flavolineata.longest.isoforms %in% c(0, 1))  %>%   
  filter(Megalopta_genalis.longest.isoforms %in% c(0, 1))  %>%          
  filter(Polistes_canadensis.longest.isoforms %in% c(0, 1))  %>%         
  filter(Polistes_dominula.longest.isoforms %in% c(0, 1))
nrow(filtered_maybe_singlecopy_df) # 3443

filter_rules <- c(0, 1, 2, 3)

# select rows if values are between 0 and 3
filtered_maybe_multiplecopy_df <- Orthogroups.GeneCount %>% filter(Apis_mellifera.longest.isoforms 
                                                                 %in% filter_rules)  %>%     
  filter(Ceratina_australensis.longest.isoforms %in% filter_rules)  %>%        
  filter(Ceratina_calcarata.longest.isoforms %in% filter_rules)  %>%            
  filter(Drosophila_melanogaster %in% filter_rules)  %>%                      
  filter(Liostenogaster_flavolineata.longest.isoforms %in% filter_rules)  %>%   
  filter(Megalopta_genalis.longest.isoforms %in% filter_rules)  %>%          
  filter(Polistes_canadensis.longest.isoforms %in% filter_rules)  %>%         
  filter(Polistes_dominula.longest.isoforms %in% filter_rules)
nrow(filtered_maybe_multiplecopy_df) # 9509

# select rows if values are between 0 and 3, present in at least 5 species
filtered_maybe_multiplecopy_df_no_outgroup <- filtered_maybe_multiplecopy_df %>%
  select(Orthogroup,
         Ceratina_australensis.longest.isoforms,
         Ceratina_calcarata.longest.isoforms,         
         Liostenogaster_flavolineata.longest.isoforms,
         Megalopta_genalis.longest.isoforms,
         Polistes_canadensis.longest.isoforms,        
         Polistes_dominula.longest.isoforms) 

# vector with species columns
species_col_vec <- c("Ceratina_australensis.longest.isoforms",
                   "Ceratina_calcarata.longest.isoforms",         
                   "Liostenogaster_flavolineata.longest.isoforms",
                   "Megalopta_genalis.longest.isoforms",
                   "Polistes_canadensis.longest.isoforms",        
                   "Polistes_dominula.longest.isoforms")

# vector will get filled by the loop with number of species without orthogene
num_absence_vec <- c()

# run through each row, obtain number of species without a orthocopy
for(i in 1:nrow(filtered_maybe_multiplecopy_df_no_outgroup)){
  num_absence_vec <- c(num_absence_vec,
                       sum(filtered_maybe_multiplecopy_df_no_outgroup[i, species_col_vec] == 0))
    
  
}

# subset for orthogroups for which a minimum of 5 species have between one and three copies
# so for when num_absence_vec <= 1: 3718 orthogroups, which is a good estimate 
# when we know that there are around 4,000 1:1 orthogroups in 6 orther species of Hymenoptera
restricted_orthogroups_df<- filtered_maybe_multiplecopy_df_no_outgroup[num_absence_vec <= 1, ]

#nrow(restricted_orthogroups_df)

# check how many copy classes there are
copy_classes_list <- apply(restricted_orthogroups_df[, species_col_vec],
      1,
      function(x) sort(unique(as.numeric(unname(x)))))

# this list needs to be filtered to show only the unqiue combinations
copy_classes_list[!duplicated(lapply(copy_classes_list, sort))]

# there are 12 classes

# no absence, only single (copy_classes_list length is 1 and sum is 1)
# 1

# no absence, only double (copy_classes_list length is 1 and sum is 2)
# 2



# absence or in single copy [length 2, sum 1]
#  0 1

# absent or in double copies [length 2, sum 2]
# 0 2

# no absence, includes single, double (copy_classes_list length is 2 and sum is 3)
# 1 2

# no absence, includes single, triple (copy_classes_list length is 2 and sum is 4)
# 1 3

# no absence, includes double, triple (copy_classes_list length is 2 and sum is 5)
# 2 3






# absent or in single or double copies [length 3, sum 3]
# 0 1 2

# absent or in single or triple copies [length 3, sum 4]
# 0 1 3

# absent or in double or triple copies [length 3, sum 5]
# 0 2 3

# no absence, includes single, double, triple (copy_classes_list length is 3 and sum is 6)
# 1 2 3


# absent or in single, double or triple copies [length 4, sum 6]
# 0 1 2 3


# need a dataframe with
# x = types of classes
# y = number of orthogroups
# colour = copy number

# how many orthogroups per class?
# make a summary table
summary_table <- as.data.frame(matrix(0, nrow = 12, ncol = 2))
colnames(summary_table) <- c("orthogroupClasses", "count")
summary_table$orthogroupClasses <- c("OnlySingle",
                                     "OnlyDouble",
                                     "AlmostSingle",
                                     "AlmostDouble",
                                     "SingleDouble",
                                     "SingleTriple",
                                     "DoubleTriple",
                                     "AlmostSingleDouble",
                                     "AlmostSingleTriple",
                                     "AlmostDoubleTriple",
                                     "SingleDoubleTrible",
                                     "AlmostSingleDoubleTrible")
# set the counter for true 1:1 orthologs at 0
OnlySingle_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 1 & sum(copy_classes_list[[i]]) == 1){
    OnlySingle_count <- OnlySingle_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "OnlySingle"] <- OnlySingle_count

#################
# set the counter for true 2:2 orthologs at 0
OnlyDouble_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 1 & sum(copy_classes_list[[i]]) == 2){
    OnlyDouble_count <- OnlyDouble_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "OnlyDouble"] <- OnlyDouble_count

#################
# set the counter for almost 1:1 orthologs (1, 0) at 0
AlmostSingle_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 2 & sum(copy_classes_list[[i]]) == 1){
    AlmostSingle_count <- AlmostSingle_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "AlmostSingle"] <- AlmostSingle_count

#################
# absent or in double copies [length 2, sum 2]
# 0 2
# set the counter for almost 2:2 orthologs (2, 0) at 0
AlmostDouble_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 2 & sum(copy_classes_list[[i]]) == 2){
    AlmostDouble_count <- AlmostDouble_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "AlmostDouble"] <-AlmostDouble_count

#################
# single or in double copies [length 2, sum 2]
# 1 2

SingleDouble_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 2 & sum(copy_classes_list[[i]]) == 3){
    SingleDouble_count <- SingleDouble_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "SingleDouble"] <-SingleDouble_count

#################
# single or in triple copies [length 2, sum 4]
# 1 2

SingleTriple_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 2 & sum(copy_classes_list[[i]]) == 4){
    SingleTriple_count <- SingleTriple_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "SingleTriple"] <-SingleTriple_count

#################
# no absence, includes double, triple (copy_classes_list length is 2 and sum is 5)
# 2 3
DoubleTriple_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 2 & sum(copy_classes_list[[i]]) == 5){
    DoubleTriple_count <- DoubleTriple_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "DoubleTriple"] <-DoubleTriple_count

#################
# absent or in single or double copies [length 3, sum 3]
# 0 1 2
AlmostSingleDouble_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 3 & sum(copy_classes_list[[i]]) == 3){
    AlmostSingleDouble_count <- AlmostSingleDouble_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "AlmostSingleDouble"] <-AlmostSingleDouble_count

#################
# absent or in single or triple copies [length 3, sum 4]
# 0 1 3
AlmostSingleTriple_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 3 & sum(copy_classes_list[[i]]) == 4){
    AlmostSingleTriple_count <- AlmostSingleTriple_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "AlmostSingleTriple"] <-AlmostSingleTriple_count

#################
# absent or in double or triple copies [length 3, sum 5]
# 0 2 3

AlmostDoubleTriple_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 3 & sum(copy_classes_list[[i]]) == 5){
    AlmostDoubleTriple_count <- AlmostDoubleTriple_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "AlmostDoubleTriple"] <-AlmostDoubleTriple_count

#################
# no absence, includes single, double, triple (copy_classes_list length is 3 and sum is 6)
# 1 2 3
SingleDoubleTrible_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 3 & sum(copy_classes_list[[i]]) == 6){
    SingleDoubleTrible_count <- SingleDoubleTrible_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "SingleDoubleTrible"] <-SingleDoubleTrible_count

#################
# absent or in single, double or triple copies [length 4, sum 6]
# 0 1 2 3
AlmostSingleDoubleTrible_count <- 0

# count number of true 1:1 orthologs
for(i in 1:length(copy_classes_list)){
  if(length(copy_classes_list[[i]]) == 4 & sum(copy_classes_list[[i]]) == 6){
    AlmostSingleDoubleTrible_count <- AlmostSingleDoubleTrible_count + 1
  }
}
summary_table$count[summary_table$orthogroupClasses == "AlmostSingleDoubleTrible"] <-AlmostSingleDoubleTrible_count


summary_table
sum(summary_table$count)

# make into a named vector
summary_expression_vec <- summary_table$count
names(summary_expression_vec) <- c("one",
                                   "two",
                                   "one&zero",
                                   "two&zero",
                                   "one&two",
                                   "one&three",            
                                   "two&three",
                                   "one&two&zero",
                                   "one&three&zero",
                                   "two&three&zero",
                                   "one&two&three",
                                   "one&two&three&zero")

upset(fromExpression(summary_expression_vec), order.by = "freq")

# save the orthogroup list
write.table(x = restricted_orthogroups_df,
      file = "2021-04-23-relaxed-filtered-orthogroup-list",
      row.names = FALSE,
      quote = FALSE)

##





#####################################
#  filter on the same rule for only bees

# select rows if values are between 0 and 3
bees_filtered_maybe_multiplecopy_df <- Orthogroups.GeneCount %>% filter(Apis_mellifera.longest.isoforms 
                                                                   %in% filter_rules)  %>%     
  filter(Ceratina_australensis.longest.isoforms %in% filter_rules)  %>%        
  filter(Ceratina_calcarata.longest.isoforms %in% filter_rules)  %>%            
  filter(Drosophila_melanogaster %in% filter_rules)  %>%                      
  filter(Megalopta_genalis.longest.isoforms %in% filter_rules)
nrow(bees_filtered_maybe_multiplecopy_df) # 9742

# select rows if values are between 0 and 3, present in at least 5 species
bees_filtered_maybe_multiplecopy_df_no_outgroup <- bees_filtered_maybe_multiplecopy_df %>%
  select(Orthogroup,
         Ceratina_australensis.longest.isoforms,
         Ceratina_calcarata.longest.isoforms,         
         Megalopta_genalis.longest.isoforms) 

# vector with species columns
bees_species_col_vec <- c("Ceratina_australensis.longest.isoforms",
                     "Ceratina_calcarata.longest.isoforms",
                     "Megalopta_genalis.longest.isoforms")

# vector will get filled by the loop with number of species without orthogene
bees_num_absence_vec <- c()

# run through each row, obtain number of species without a orthocopy
for(i in 1:nrow(bees_filtered_maybe_multiplecopy_df_no_outgroup)){
  bees_num_absence_vec <- c(bees_num_absence_vec,
                       sum(bees_filtered_maybe_multiplecopy_df_no_outgroup[i, bees_species_col_vec] == 0))
  
  
}

# subset for orthogroups for which a minimum of 2 species have between one and three copies
# so for when num_absence_vec <= 1: 5787 orthogroups
bees_restricted_orthogroups_df<- bees_filtered_maybe_multiplecopy_df_no_outgroup[bees_num_absence_vec <= 1, ]


# save the orthogroup list
write.table(x = bees_restricted_orthogroups_df,
            file = "2021-04-26-relaxed-filtered-orthogroup-list-bees",
            row.names = FALSE,
            quote = FALSE)


##################################
# just for wasps
# select rows if values are between 0 and 3
wasps_filtered_maybe_multiplecopy_df <- Orthogroups.GeneCount %>% filter(Apis_mellifera.longest.isoforms 
                                                                   %in% filter_rules)  %>%     
  filter(Drosophila_melanogaster %in% filter_rules)  %>%                      
  filter(Liostenogaster_flavolineata.longest.isoforms %in% filter_rules)  %>%   
  filter(Polistes_canadensis.longest.isoforms %in% filter_rules)  %>%         
  filter(Polistes_dominula.longest.isoforms %in% filter_rules)
nrow(wasps_filtered_maybe_multiplecopy_df) # 10083

# select rows if values are between 0 and 3, present in at least 5 species
wasps_filtered_maybe_multiplecopy_df_no_outgroup <- wasps_filtered_maybe_multiplecopy_df %>%
  select(Orthogroup,
         Liostenogaster_flavolineata.longest.isoforms,
         Polistes_canadensis.longest.isoforms,        
         Polistes_dominula.longest.isoforms) 

# vector with species columns
wasps_species_col_vec <- c("Liostenogaster_flavolineata.longest.isoforms",
                     "Polistes_canadensis.longest.isoforms",        
                     "Polistes_dominula.longest.isoforms")

# vector will get filled by the loop with number of species without orthogene
wasps_num_absence_vec <- c()

# run through each row, obtain number of species without a orthocopy
for(i in 1:nrow(wasps_filtered_maybe_multiplecopy_df_no_outgroup)){
  wasps_num_absence_vec <- c(wasps_num_absence_vec,
                       sum(wasps_filtered_maybe_multiplecopy_df_no_outgroup[i,
                                                                            wasps_species_col_vec] == 0))
  
  
}

# subset for orthogroups for which a minimum of 2 species have between one and three copies
# so for when num_absence_vec <= 1: 6983
wasps_restricted_orthogroups_df<- wasps_filtered_maybe_multiplecopy_df_no_outgroup[wasps_num_absence_vec <= 1, ]

# save the orthogroup list
write.table(x = wasps_restricted_orthogroups_df,
            file = "2021-04-26-relaxed-filtered-orthogroup-list-wasps",
            row.names = FALSE,
            quote = FALSE)
