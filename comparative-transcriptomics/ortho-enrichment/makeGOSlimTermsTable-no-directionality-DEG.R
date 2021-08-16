#################################
#### GO Slim pipeline Step 2 ####
#################################

# load  libraries
library("biomaRt")
library("topGO")

##########################
# vector of all species available
species_vec <- c("Ceratina_australensis",
                 "Ceratina_calcarata",
                 "Liostenogaster_flavolineata",
                 "Megalopta_genalis",
                 "Polistes_canadensis",
                 "Polistes_dominula")

# vector of all categories available
goCategory_vec <- c("BP", "MF", "CC")

# loop over the species
for(myspecies in species_vec){
  # loop over the categories
  for(gocat in goCategory_vec){
    
    # set the arguments
    args <- c(myspecies, 
              gocat)
    
    # set experiment details (species, GO Category)
    this_species     <- args[1]
    this_goCategory  <- args[2]

    # set species data input
    raw_results_file         <- paste("inputOnScratch/", this_species,
                                      "_DESeq2_results.txt", sep = "")
    raw_selectionGenes_file  <- paste("inputOnScratch/", this_species,
                                      "_all_DEGs.txt", sep = "")
    raw_blast_results_file   <- paste("resultOnScratch/", this_species,
                                      "_filtered", sep = "")
    hash_table <- paste("tmp/", this_species,
                        "_protein_gene_hash_table", sep = "")
    
    goslimmer_result_file <- paste("inputOnScratch/",
                                   this_species,
                                   "-FlybaseGeneID-GOSlimID.tabular", sep = "")
    # import all DEG results 
    # columns: X baseMean log2FoldChange lfcSE stat pvalue padj
    # rows: genes
    raw_results <- read.delim(raw_results_file,
                              stringsAsFactors = FALSE)
    
    # import significant DEG 
    # one column with names of genes
    raw_selectionGenes <- read.table(raw_selectionGenes_file,
                                     quote = "\"",
                                     comment.char = "",
                                     stringsAsFactors = FALSE)
    
    
    # blast results: each Drosophila Orthogroup has a match in similarity with a species' protein
    raw_blast_results <- read.table(raw_blast_results_file,
                                    stringsAsFactors = FALSE)
    
    # hash table: protein in column 1, gene in column 2
    hash_df <- read.table(hash_table,
                          stringsAsFactors = FALSE)
    
    # GO Slim Term results
    # flybase gene id and GOSlim Terms
    # genes whose protein sequences match C. australensis protein
    this_species.FlybaseGeneID.GOSlimID <- read.delim(goslimmer_result_file,
                                                   stringsAsFactors = FALSE)
    
    
    # add column names 
    colnames(raw_results) <- c("gene", colnames(raw_results)[2:7])
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
    
    # table with translation name, gene name, full set GO ID
    gene2Go_raw3 <- getBM(attributes = c("flybase_translation_id",
                                         "flybase_gene_id",
                                         "go_id"), 
                          filters     = "flybase_translation_id", 
                          values      = droso_gene_list, 
                          mart        = ensembl,
                          useCache    = FALSE)
    
    # make a list with flybase translation id (element) and GO slim terms (vector in each element)
    this_species.FlybaseTranslationID.GOSlimID <- vector(mode = "list",
                                                                  length = nrow(gene2Go_raw3))
    
    # add named elements (translation ID)
    names(this_species.FlybaseTranslationID.GOSlimID) <- gene2Go_raw3$flybase_translation_id
    
    # make the equivalent of gene2go (the annotation input to TopGO new)
    # fill each element with GOSlim terms
    
    for(i in 1:length(names(this_species.FlybaseTranslationID.GOSlimID))){
      # loop through each flybase transcript id e.g. "FBpp0070031"
      #i <- 1
      
      # gene id corresponding to the translation id e.g. "FBgn0031086"
      gene_id <- gene2Go_raw3$flybase_gene_id[match(names(this_species.FlybaseTranslationID.GOSlimID)[i],
                                                    gene2Go_raw3$flybase_translation_id)]
      
      # GOSlim terms corresponding to the gene_id e.g. "GO:0030154"
      GOSlim_term_vec <- this_species.FlybaseGeneID.GOSlimID$GOSlim.Term[this_species.FlybaseGeneID.GOSlimID$Gene == gene_id]
      
      # add it to the table
      this_species.FlybaseTranslationID.GOSlimID[[i]] <- GOSlim_term_vec
    }
    
    # this data is now ready for TopGo
    
    
    ##### GO term enrichment 
    
    # aim to change species's protein names for drosophila names
    # because the TopGO database does not contain non-model data
    
    # there are NA because blasting droso against the species might have produced no hit
    raw_results$droso_gene <- raw_blast_results$sseqid[match(raw_results$gene,
                                                             raw_blast_results$qseqid)]
    
    raw_selectionGenes$droso_gene <- raw_blast_results$sseqid[match(raw_selectionGenes$gene,
                                                                    raw_blast_results$qseqid)]
    
    
    # remove NA. 6,291 species gene have an droso hit, including 247 in the selection
    raw_results <- raw_results[!is.na(raw_results$droso_gene), ]
    raw_selectionGenes <- raw_selectionGenes[!is.na(raw_selectionGenes$droso_gene), ]
    
    
    ## make a vector with 0 or 1 values depending if a gene is DE or not
    # results: lists of genes differentially expressed 
    geneList <- rep(0, times = length(rownames(raw_results)))
    
    # name each value with the droso genes names
    names(geneList) <- raw_results$droso_gene
    
    # selectionGenes: list of DEG for selection
    DEGenes <- raw_selectionGenes$droso_gene
    
    # for each gene that is the focus of the analysis, change the value 0 for 1
    geneList[DEGenes] <- 1
    
    # change the class to factor
    geneList <-  as.factor(geneList)
    
    
    ## Build the topGO object for biological process ontology
    this_topGOdata <- new("topGOdata",
                          ontology = this_goCategory,
                          allGenes = geneList,
                          geneSel  = DEGenes,
                          nodeSize = 5,
                          annot    = annFUN.gene2GO,
                          gene2GO  = this_species.FlybaseTranslationID.GOSlimID)
    
    # test for enrichment
    # because we coded the genes 1 or 0 for DEG presence or absence, Fisher test (gene count) is probably the best algorithm
    # classic: each GO category is tested independently
    this_topGOresult <- runTest(this_topGOdata,
                                algorithm = "classic",
                                statistic = "fisher")
    
    
    
    # create a result table
    # GO Terms identified by fisher test
    myTable <- GenTable(this_topGOdata,
                        pvalue = this_topGOresult,
                        topNodes = length(this_topGOdata@graph@nodes),
                        numChar = 100)
    
    
    # add columns to specify test details
    myTable$species <- this_species
    myTable$goCategory <- this_goCategory
    
    # make a file name
    this_file_name <- paste("resultOnScratch/topgo_result_slimTerms",
                            this_species,
                            this_goCategory,
                            sep = "_")
    # save table
    write.table(x          = myTable,
                file       = this_file_name,
                quote      = FALSE,
                row.names  = FALSE,
                sep        = "\t")
    
    
  }
}
    
    
