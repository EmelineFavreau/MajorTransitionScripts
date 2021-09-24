#################################
#### GO Slim pipeline Step 1 ####
#################################

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
  raw_blast_results_file   <- paste("resultOnScratch/", this_species,
                                    "_filtered", sep = "")
  hash_table <- paste("tmp/", this_species,
                      "_protein_gene_hash_table", sep = "")
  
  
  
  # blast results: each Drosophila Orthogroup has a match in similarity with a species' protein
  raw_blast_results <- read.table(raw_blast_results_file,
                                  stringsAsFactors = FALSE)
  
  # hash table: protein in column 1, gene in column 2
  hash_df <- read.table(hash_table,
                        stringsAsFactors = FALSE)
  
  # add column names 
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
  # for go slim I need the flybase gene id, not the protein id
  gene2GoSlim_raw <- getBM(attributes = c("flybase_gene_id", "go_id"), 
                        filters     = "flybase_translation_id", 
                        values      = droso_gene_list, 
                        mart        = ensembl,
                        useCache    = FALSE)
  
  # Remove the genes without GO terms
  gene2GoSlim_df <- subset(x = gene2GoSlim_raw,
                        subset = !go_id == "")
  # vec of file name
  fileName <- paste(myspecies, "gene2GoSlim.txt", sep = "_")
  
  # save it to import into Galaxy GoSlimmer
  write.table(gene2GoSlim_df,
              file      = fileName,
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE,
              sep       = "\t")
}
### END OF SCRIPT ONE


