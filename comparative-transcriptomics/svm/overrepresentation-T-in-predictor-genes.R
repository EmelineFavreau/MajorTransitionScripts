
# https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
# we have 4145 genes, 114 are predictor
# we have 10 predictor genes that are TF, 
total_num_genes <- 4145
total_predictor_genes <- 114
TF_in_predictor_gene_num <- 10
TF_in_non_predictor_gene_num <- 628
predictor_gene_num <- total_predictor_genes - TF_in_predictor_gene_num
non_predictor_gene_num <- total_num_genes - TF_in_non_predictor_gene_num

# focus on total TF (not focused on evidence-based)
# Is there an overrepresentation of TF in predictor genes than in all orthogroups?
# make a contigency table
d <- data.frame(
                In_predictor_genes = c(TF_in_predictor_gene_num,
                                       predictor_gene_num),
                not_in_predictor_genes = c(TF_in_non_predictor_gene_num,
                                           non_predictor_gene_num))

row.names(d) <- c("TF", "non_TF")

# contigency table
d

# odd ratio is how many more TF in predictor genes than in non-predictor genes
# 0.3111847 


# one-sided test 
# p1 is the first column
# p2 is the second column
# comparing the null hypothesis p1=p2
# alternative hypothesis p1>p2 (overrepresentation of TF in predictor genes)
fisher.test(d, alternative = "greater")

#p-value = 0.9836
# no over-representation of transcription factors in the predictor genes