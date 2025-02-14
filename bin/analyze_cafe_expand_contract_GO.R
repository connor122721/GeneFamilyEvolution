# Read significant trees from CAFE5 and do GO term enrichment testing
# Connor Murray 1.23.2024; modified 2.1.2025
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(AnnotationDbi)
library(phytools)
library(ape)
library(tidyverse)
library(foreach)
library(ggtree)
library(viridis)
library(clusterProfiler)
library(parallel)
library(argparse)
cores <- 3
cl <- parallel::makeCluster(cores)

# Common theme element for plotting
themei <- (
  theme_bw() +
    theme(strip.text = element_text(face="bold", size=18),
          legend.text = element_text(face="bold", size=18),
          title = element_text(face="bold", size=20),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20)))

# Add argparse to parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--cafe_dir", default="output/cafe", help="Path to CAFE result directory")
args <- parser$parse_args()

# Use the parsed argument instead of the hardcoded value
cafe_dir <- args$cafe_dir
# cafe_dir = "/project/berglandlab/connor/GeneFamilyEvolution/output/cafe/"

# Define directory and list CAFÉ result text files
result_files <- list.files(cafe_dir, "results.txt$", recursive = TRUE, full.names = TRUE)

# Function to extract -lnL from a file
extract_lnL <- function(file) {
  lines <- readLines(file)
  lnL_line <- grep("-lnL", lines, value = TRUE)
  if (!length(lnL_line)) lnL_line <- grep("lnL", lines, value = TRUE)
  if (!length(lnL_line)) return(NA_real_)
  as.numeric(sub(".*?(-?\\d+\\.?\\d*).*", "\\1", lnL_line[1]))
}

# Extract likelihoods and identify the best run
lnL_vals <- sapply(result_files, extract_lnL)
if (all(is.na(lnL_vals))) stop("No -lnL values found in any .txt file.")
best_file <- result_files[which.min(lnL_vals)]
cat("Best run:", best_file, "with -lnL =", lnL_vals[which.min(lnL_vals)], "\n")

# Load significant trees from the best run folder
tre_file <- file.path(dirname(best_file), "Significant_trees.tre")
if (!file.exists(tre_file)) stop("File not found:", tre_file)
tre <- read.nexus(tre_file)

# Hog -> OG metadata
og <- data.table(fread("/project/berglandlab/connor/GeneFamilyEvolution/output/longest_orf/primary_transcripts/OrthoFinder/Results_Feb04/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"))
colnames(og)[3] <- "Gene_Tree" 

# Hogs annotation
hog.ann <- data.table(fread("/project/berglandlab/connor/GeneFamilyEvolution/output/orthofinder/hogs.function.tsv") %>% 
                        left_join(og, by=c("Gene"="HOG")))
colnames(hog.ann)[4] <- "Gene_Tree" 

# A helper to extract the gene count from a CAFE5 annotated tree.
extract_count <- function(label) {
  as.numeric(gsub(".*_", "", label))
}

# New annotate_tree function: 
annotate_tree <- function(tree.i) {
  
  # tree.i is assumed to be a "phylo" object from CAFE5.
  ntip <- length(tree.i$tip.label)
  
  # Build a data.table with the tip labels and tip indices.
  tip_dt <- data.table(File = tree.i$tip.label,
                       tip_index = 1:ntip)
  
  # Extract species name and gene count from the tip label.
  # For example, a tip like "anopheles<1>_2" will yield spp = "anopheles" and tip_count = 2.
  tip_dt <- data.table(tip_dt %>%
    mutate(spp = tstrsplit(File, "<", fixed = TRUE)[[1]],
           tip_count = as.numeric(tstrsplit(File, "_", fixed = TRUE)[[2]])))
  
  # For each tip, find the parent node from the tree's edge matrix.
  # The edge matrix has two columns: [parent, child].
  # We match the tip index (child) and extract the parent's index.
  tip_dt[, parent_index := tree.i$edge[match(tip_index, tree.i$edge[,2]), 1]]
  
  # Get the parent's label.
  # In ape objects, tip indices are 1...ntip and node indices are (ntip+1):(ntip+nnode).
  # If parent's index > ntip, then parent's label is in tree.i$node.label at position (parent_index - ntip).
  tip_dt[, parent_label := ifelse(parent_index > ntip,
                                  tree.i$node.label[parent_index - ntip],
                                  tree.i$tip.label[parent_index])]
  
  # Extract the parent's gene count.
  tip_dt[, parent_count := extract_count(parent_label)]
  
  # Compute the difference: (tip gene count) minus (parent gene count).
  tip_dt[, change_delt := tip_count - parent_count]
  
  # Label the change: "E" for expansion (if tip count is higher), "C" for contraction (if lower),
  # and "S" for same.
  tip_dt[, change_delt_col := ifelse(change_delt > 0, "E",
                                     ifelse(change_delt < 0, "C", "S"))]
  
  # Optionally, mark if the original tip label contains a "*" (as in your previous code).
  tip_dt[, sig := ifelse(grepl("\\*", File), "*", "")]
  
  return(tip_dt)
}

# Extract name
tree.ann <- mclapply(tre, annotate_tree, mc.cores = cores)
tree.dt <- data.table(rbindlist(tree.ann, use.names = T, fill=T, idcol = "tree"))
tree.tally <- data.table(tree.dt %>%
                group_by(tree) %>%
                count(sig) %>% 
              mutate(SigNum=case_when(sig==""~"N",
                                      sig=="*"~"Y")))

# Count NUmber of significant changes and add metadata
tree.dt1 <- data.table(tree.dt %>% 
              group_by(tree) %>%
              summarize(mean=mean(change_delt)) %>%
              left_join(tree.tally)  %>% 
              left_join(hog.ann, by=c("tree"="Gene")))

#### EXPANDED AND CONTRACTED GENES ####

# Expanded gene families
expand <- data.table(tree.dt[change_delt>0][sig=="*"] %>%
                       left_join(tree.tally)  %>% 
                       left_join(hog.ann, by=c("tree"="Gene")))

# Extract expanding genes
extract.expanding.genes <- function(spp.i) {
  #spp.i="magna"; 
  data=expand
  spp.out <- paste(spp.i,"protein",sep=".")
  ppp <- select(data[spp==spp.i], contains(spp.out)) %>% 
    drop_na() %>%
    separate_rows({{  spp.out  }}, sep = ", ") %>%
    unlist() %>% unname() %>% unique()
  
  write.table(ppp, file = paste(spp.out, ".expanding.genes.txt", sep=""), 
              quote = F, row.names = F, col.names = F)
  
  return(data.table(ppp, species=spp.i))
}

# Find expanding gene families
ext.list <- lapply(unique(tree.dt$spp), FUN=extract.expanding.genes)

# Contracted genes
contract <- data.table(tree.dt[change_delt<0][sig=="*"] %>%
                         left_join(tree.tally)  %>% 
                         left_join(hog.ann, by=c("tree"="Gene")))

# Extract contracting genes
extract.contract.genes <- function(spp.i) {
  #spp.i="magna"; 
  data=contract
  spp.out <- paste(spp.i,"protein",sep=".")
  ppp <- select(data[spp==spp.i], contains(spp.out)) %>% 
    drop_na() %>%
    separate_rows({{  spp.out  }}, sep = ", ") %>%
    unlist() %>% unname() %>% unique()
  
  write.table(ppp, file = paste(spp.out, ".contract.genes.txt", sep=""), 
              quote = F, row.names = F, col.names = F)
  
  return(data.table(ppp, species=spp.i))
}

con.list <- lapply(unique(tree.dt$spp), FUN=extract.contract.genes)

# Save output
saveRDS(contract, file = "contracted.genes.rds")
saveRDS(expand, file = "expanded.genes.rds")

############## ENRICHMENT ###############

# read in table of genes belonging to orthogroups, select significant OGs and parse into gene lists
gene_lists_tbl_contract <- og %>%
  filter(HOG %in% contract[Gene_Tree=="n0"]$tree) %>%
  dplyr::select(-c(HOG,OG,Gene_Tree))

gene_lists_tbl_expand <- og %>%
  filter(HOG %in% expand[Gene_Tree=="n0"]$tree) %>%
  dplyr::select(-c(HOG,OG,Gene_Tree))

# function to extract gene list from gene_lists_tbl
get_gene_list <- function(genome = "string"){
  gene_lists_tbl[genome] %>% 
    drop_na() %>%
    separate_rows({{  genome  }}, sep = ", ") %>%
    unlist() %>% unname() %>% unique()
}

# Annotation files from Blast2GO/OmicsBox
GOtables_filenames <- list.files(path = "/project/berglandlab/connor/GeneFamilyEvolution/output/orthogroups_annotation", 
                                 pattern = "GOterm_mapping.tsv",  
                                 full.names = TRUE)

# Parsing function (file to table)
GOtable_parsing <- function(infile = blast2go.tsv) {
  # Extract Genome name from file name (specific to this file name structure)
  #infile=GOtables_filenames[3]
  
  GOtable_genome <- infile %>% 
    str_remove(".*/") %>%
    str_remove("_.*")
  
  # read table from file and parse
  GOtable <- data.table(read.delim(infile, header=F)) %>%
    mutate(Genome = GOtable_genome) %>%
    dplyr::rename(GeneID = "V1") %>%
    dplyr::rename(GO_ID = "V2") %>%
    separate_rows(GO_ID, sep = ",") %>%
    mutate(GO_ID = str_trim(GO_ID),
           GO_ID = str_replace(GO_ID, ".:GO:", "GO:")) %>%
    unite("GeneID", Genome, GeneID, sep = "_", remove = FALSE)
  
  return(GOtable)
}

# Create a simplified GO data frame
GeneSetCollection_constructor <- function(genome_GO_table = GOtable) {
  #genome_GO_table=GO_tables[[1]]
  
  go_df_simple <- genome_GO_table %>%
    dplyr::select(gene = GeneID, GO = GO_ID) %>%
    filter(!is.na(GO)) %>%
    as.data.frame()
  
  # Extract species name from the Genome column
  go_species <- genome_GO_table %>%
    dplyr::select(Genome) %>% 
    unique() %>% 
    unlist() %>% 
    unname()
  
  # Build and save the GO map
  gmap <- buildGOmap(go_df_simple)
  saveRDS(gmap, file = paste0("gmap.", go_species, ".rds"))
  
  # Create the full gene universe.
  full_universe <- genome_GO_table %>%
    dplyr::select(GeneID) %>%
    unique() %>%
    unlist() %>%
    unname()
  
  print(go_species)
  
  return(list(Species = go_species, 
              Universe = full_universe))
}

# Parse GO data table and make ready for clusterprofiler
GO_tables <- mclapply(GOtables_filenames, GOtable_parsing, mc.cores = cores)
full_gene_sets <- mclapply(GO_tables, GeneSetCollection_constructor, mc.cores = cores)
names(full_gene_sets) <- purrr::map_chr(full_gene_sets, ~ .x$Species[1])

print("Running GOTerm Enrichment!")

# Enrichment analysis function
go.test.clusterpro <- function(full_set, gene_lists_tbl, namei) {
  # full_set=full_gene_sets[[6]]; gene_lists_tbl=gene_lists_tbl_expand[[6]]; namei="expand"
  
  # Extract species ID (assumes full_set contains a Species field)
  go_species <- data.table(gene = full_set$Universe) %>%
    mutate(spp = full_set$Species) %>% 
    dplyr::select(spp) %>% unique() %>% unlist() %>% unname()
  
  print(go_species)
  print(namei)
  
  # Build candidate gene list using the appropriate column from gene_lists_tbl
  gene.cand <- gene_lists_tbl %>% 
    unlist() %>% 
    strsplit(split = ",") %>%
    unlist() %>% 
    paste0(go_species, "_", .) %>% 
    gsub(pattern = " ", replacement = "") %>% 
    unique()
  
  # Universe: the set of longest transcript proteins (unique)
  univ <- unique(full_set$Universe)
  
  # Load species-specific gene-to-GO mapping file
  spp_gmap_file <- list.files(pattern = paste0("gmap.", go_species, ".rds"), full.names = TRUE)
  if (length(spp_gmap_file) == 0) {
    stop("No gene mapping file found for species: ", go_species)}
  
  spp.gmap <- readRDS(spp_gmap_file[1])
  
  # Enrichment test using clusterProfiler's enricher()
  go.test <- enricher(gene = gene.cand,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.01,
                      universe = univ, 
                      TERM2GENE = spp.gmap %>% dplyr::select(GO, gene))
  
  go.test2 <- data.table(go.test@result) %>% filter(p.adjust < 0.05)
  
  # Iterate over each enriched GO term to compute odds ratios and add term descriptions
  pp <- foreach(i = seq_len(nrow(go.test2)), .combine = "rbind", .errorhandling = "remove") %do% {
    current_id <- go.test2$ID[i]
    message(paste(i, current_id, sep = " | "))
    go.test2[ID == current_id] %>% 
      mutate(term = go2term(ID)$Term,          # Make sure go2term() is defined
             x = eval(parse(text = GeneRatio)),
             y = eval(parse(text = BgRatio)),
             odd = x / y)}
  
  # Plot enriched GO terms
  plot.go <- { pp %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(p.adjust)) %>% 
    ggplot(aes(x = log2(odd),
               y = reorder(term, log2(odd)),
               size = Count,
               color = log10(p.adjust))) +
    geom_point() +
    labs(x = "log2(Enrichment)", 
         y = "",
         color = "FDR P-value",
         size = "Candidate gene number") +
    theme_bw() +
    scale_color_viridis(option = "magma") +
    theme(title = element_text(face = "bold", size = 15),
          legend.text = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 16),
          legend.background = element_blank(),
          axis.text.x = element_text(face = "bold", size = 15),
          axis.text.y = element_text(face = "bold", size = 15),
          axis.title.x = element_text(face = "bold", size = 18),
          axis.title.y = element_text(face = "bold", size = 18))
  }
  
  ggsave(plot = plot.go, 
    filename = paste0(namei, ".", go_species, ".goterm.cafe.pdf"), 
    width = 20, height = 12)
  
  # Write output for REVIGO
  fwrite(pp %>% select(ID, p.adjust, odd), 
         file = paste0(go_species, ".", namei, ".revigo.cafe.txt"), 
         quote = FALSE, sep = "\t")
  
  return(pp)
}

go.out.con <- lapply(full_gene_sets,
                     go.test.clusterpro,
                     gene_lists_tbl = gene_lists_tbl_contract,
                     namei = "contract")

go.out.exp <- lapply(full_gene_sets,
                     go.test.clusterpro,
                     gene_lists_tbl = gene_lists_tbl_expand,
                     namei = "expanded")

# Add GO term reduction
revigo_list_contract <- list.files(pattern = "contract.revigo")
revigo_list_expand <- list.files(pattern = "expand.revigo")

# Semantic reduction
library(rrvgo)

# Function for redundancy
reduce.go <- function(go_list_enriched) {
  # gene_list_enriched=revigo_list_contract[[1]]
  
  # Extract species ID
  geneii=data.table(fread(gene_list_enriched))
  colnames(geneii) <- c("ID", "qvalue", "OR")
  go_species <- tstrsplit(gene_list_enriched, ".", fixed=T)[[1]]
  print(go_species)
  
  # Conditional based on length
  if(nrow(geneii)<15) {
    print("List is very short, no redundancy analysis preformed.")
    reduced.data <- data.table(geneii, species=go_species) 
  } else {
    
    # Sim reduction
    simMatrix <- calculateSimMatrix(geneii$ID,
                                    orgdb="org.Dm.eg.db")
    
    # Fin reduced
    scores <- setNames(-log10(geneii$qvalue), 
                       geneii$ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.Dm.eg.db")
    
    # Collapse data
    reduced.data <- data.table(reducedTerms %>% 
                                 left_join(geneii),
                               species=go_species)
    
  }
  
  # Finish function
  return(reduced.data)
  
}

# Output redueced GO terms - reduce based on Drosophila melanogaster GO mapping
#go.reduce.out.con <- mclapply(go.out.con, reduce.go, mc.cores = cores)

# Finalize output
#go.reduce.fin.con <- data.table(rbindlist(go.reduce.out.con, fill=T))

# Save output
#saveRDS(go.reduce.out.con, file = "reduced.terms.contracting.rds")
#saveRDS(go.reduce.out.con, file = "reduced.terms.expanding.rds")