#########################################
#### GO Slim PANTHER pipeline Step 1 ####
#########################################

# specifically, we need the Dmel genes id associated with WGCNA proteins
# and the Dmel gene id associated with all proteins of given species (blasted)

# load  libraries
library("biomaRt")
library("topGO")

##########################
# 1 - input species data, blast result,
# query biomart to obtain drosophila translation names 


# vector of all species available
species_vec <- c("Ceratina_australensis",
                 "Ceratina_calcarata",
                 "Liostenogaster_flavolineata",
                 "Megalopta_genalis",
                 "Polistes_canadensis",
                 "Polistes_dominula")



#myspecies <- "Ceratina_australensis"
# loop over the species
for(myspecies in species_vec){  
  
  # set experiment details (species)
  this_species     <- myspecies
  
  # set species data input
  raw_results_file         <- paste("inputOnScratch/", this_species,
                                    "_WGCNA_Input_Genes.txt", sep = "")
  raw_selectionGenes_file  <- paste("inputOnScratch/", this_species,
                                    "_WGCNA_Significant_Module_Genes.txt", sep = "")
  raw_blast_results_file   <- paste("resultOnScratch/", this_species,
                                    "_filtered", sep = "")
  hash_table <- paste("tmp/", this_species,
                      "_protein_gene_hash_table", sep = "")
  
  
  # import WGCNA results 
  # column: genes
  raw_results <- read.delim(raw_results_file,
                            stringsAsFactors = FALSE)
  
  # import genes in significant modules (correlated to reproductive status) 
  # one column with names of genes
  raw_selectionGenes <- read.table(raw_selectionGenes_file,
                                   quote = "\"",
                                   comment.char = "",
                                   stringsAsFactors = FALSE)
  
  # blast results: each Drosophila Orthogroup has a match in similarity with
  # a species' protein
  raw_blast_results <- read.table(raw_blast_results_file,
                                  stringsAsFactors = FALSE)
  
  # hash table: protein in column 1, gene in column 2
  hash_df <- read.table(hash_table,
                        stringsAsFactors = FALSE)
  
  # add column names 
  colnames(raw_results) <- c("gene")
  colnames(raw_selectionGenes) <- "gene"
  colnames(raw_blast_results) <- c("qseqid", "sseqid", "pident", "length",
                                   "mismatch", "gapopen", "qstart",
                                   "qend", "sstart", "send", "evalue",
                                   "bitscore")
  
  colnames(hash_df) <- c("qseqid", "gene")
  
  # update blast query sequence id (to gene-LOCXXX, matching DESeq2 result table)
  raw_blast_results$qseqid <- hash_df$gene[match(raw_blast_results$qseqid,
                                                 hash_df$qseqid)]
  
  # r getting gene to go mapping droso
  
  # connect to the genes services
  ensembl <- useEnsembl(biomart = "ensembl",
                        dataset = "dmelanogaster_gene_ensembl")

  
  # list of droso genes that I want, I think these are transcript id
  droso_gene_list <- raw_blast_results$sseqid
  
  
  ##########################
  # 2 - obtain annotation file for Galaxy input (R)
  # GOSlimmer needs three files
  # full Gene Ontology file (from website, made by flybase)
  # slim gene ontology file (from website, made by flybase)
  # a tabular annotation file in two-column with gene product ids in first column, GO Terms in second one
  # for go slim pnather I need the flybase gene id, and the protein id
  gene2GoSlim_raw <- getBM(attributes = c("flybase_translation_id", 
                                          "flybase_gene_id", 
                                          "go_id"), 
                           filters     = "flybase_translation_id", 
                           values      = droso_gene_list, 
                           mart        = ensembl,
                           useCache    = FALSE)
  
  # obtain Dmel protein that are WGCNA in species
  Dmel_protein_WGCNA <- raw_blast_results$sseqid[match(raw_selectionGenes$gene,
                                                      raw_blast_results$qseqid)]
  Dmel_protein_WGCNA <- Dmel_protein_WGCNA[!is.na(Dmel_protein_WGCNA)]
  
  # obtain Dmel genes that are WGCNA in species
  Dmel_genes_WGCNA <- gene2GoSlim_raw$flybase_gene_id[match(Dmel_protein_WGCNA,
                                        gene2GoSlim_raw$flybase_translation_id)]
  
  # vec of file name
  fileName <- paste(myspecies, "Dmel_genes_WGCNA_SlimPanther.txt", sep = "_")
  
  # save it to import into Galaxy GoSlimmer
  write.table(Dmel_genes_WGCNA,
              file      = fileName,
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE,
              sep       = "\t")
  
  
  # obtain Dmel genes that are in species
  Dmel_protein <- raw_blast_results$sseqid[match(raw_results$gene,
                                                      raw_blast_results$qseqid)]
  Dmel_protein <- Dmel_protein[!is.na(Dmel_protein)]
  
  # obtain Dmel genes that are in species
  Dmel_genes <- gene2GoSlim_raw$flybase_gene_id[match(Dmel_protein,
                                                           gene2GoSlim_raw$flybase_translation_id)]
  
  # vec of file name
  fileName <- paste(myspecies, "Dmel_genes_SlimPanther.txt", sep = "_")
  
  # save it to import into Galaxy GoSlimmer
  write.table(Dmel_genes,
              file      = fileName,
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE,
              sep       = "\t")
}
### END OF SCRIPT ONE


