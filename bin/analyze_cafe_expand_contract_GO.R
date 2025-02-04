# Read significant trees from CAFE and do GO term enrichment testing
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

# Read significant trees
setwd("/project/berglandlab/connor/GeneFamilyEvolution/")
tre <- read.nexus(file="output/cafe/r1_hogs/Significant_trees.tre")

# Hog -> OG metadata
og <- data.table(fread("output/longest_orf/primary_transcripts/OrthoFinder/Results_Feb01/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"))

# Hogs annotation
hog.ann <- data.table(fread("output/orthogroups_annotation/hogs.function.tsv") %>% 
                        left_join(og, by=c("Gene"="HOG")))

# Open tree data and add change
annotate_tree <- function(tree.i) {
  # tree.i=tre[[12]]

  # Read in tree
  tree <- data.table(File=as.character(tree.i$tip.label))
  
  # Extract naming
  tree1 <- data.table(tree %>%
                     mutate(spp = tstrsplit(File,"<", fixed=T)[[1]],
                            num = tstrsplit(File, "_", fixed=T)[[2]]))
  
  tree2 <- tree1 %>% 
    mutate(sig=case_when(File %in% tree1$File[grep(tree1$File, pattern="*", fixed=T)] ~ "*",
                         TRUE ~ ""))
  
  nodei=abs(as.numeric(gsub(x = tree.i$node.label, pattern = "*", "", fixed=T)))
  
  tree3 <- tree2 %>% 
    mutate(change.delt = case_when(spp %in% c("pulicaria", "pulexnam")~as.numeric(num)-nodei[5],
                                   spp %in% c("sinensis", "magna")~as.numeric(num)-nodei[6],
                                   spp %in% c("galeata")~as.numeric(num)-nodei[3],
                                   spp %in% c("pulexeuro")~as.numeric(num)-nodei[4],
                                   spp %in% c("carinata")~as.numeric(num)-nodei[1])) %>% 
    mutate(change.delt.col=case_when(spp %in% c("pulicaria", "pulexnam") & as.numeric(num) > nodei[5] ~ "E",
                                     spp %in% c("pulicaria", "pulexnam") & as.numeric(num) < nodei[5] ~ "C",
                                     spp %in% c("sinensis", "magna") & as.numeric(num) > nodei[6] ~ "E",
                                     spp %in% c("sinensis", "magna") & as.numeric(num) < nodei[6] ~ "C",
                                     spp %in% c("galeata") & as.numeric(num) > nodei[3] ~ "E",
                                     spp %in% c("galeata") & as.numeric(num) < nodei[3] ~ "C",
                                     spp %in% c("pulexeuro")& as.numeric(num) > nodei[4] ~ "E",
                                     spp %in% c("pulexeuro") & as.numeric(num) < nodei[4] ~ "C",
                                     spp %in% c("carinata") & as.numeric(num) > nodei[1] ~ "E",
                                     spp %in% c("carinata") & as.numeric(num) > nodei[1] ~ "C",
                                     TRUE ~ "S"))
  
}

# Extract name
tree.ann <- lapply(tre, annotate_tree)
tree.dt <- data.table(rbindlist(tree.ann, use.names = T, fill=T, idcol = "tree"))
tree.tally <- data.table(tree.dt %>%
                group_by(tree) %>%
                count(sig) %>% 
              mutate(SigNum=case_when(sig==""~"N",
                                      sig=="*"~"Y")))

# Count NUmber of significant changes and add metadata
tree.dt1 <- data.table(tree.dt %>% 
              group_by(tree) %>%
              summarize(mean=mean(change.delt)) %>%
              left_join(tree.tally)  %>% 
              left_join(hog.ann, by=c("tree"="Gene")))

#### EXPANDED AND CONTRACTED GENES ####

# Expanded gene families
expand <- data.table(tree.dt[change.delt>0][sig=="*"] %>%
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
contract <- data.table(tree.dt[change.delt<0][sig=="*"] %>%
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


# saveRDS(contract, file = "/scratch/csm6hg/data/contracted.genes.sep7.rds")
# saveRDS(expand, file = "/scratch/csm6hg/data/expanded.genes.sep7.rds")

# Read in data
# contract <- data.table(readRDS("contracted.genes.sep7.rds"))
# expand <- data.table(readRDS("expanded.genes.sep7.rds"))

############## ENRICHMENT ###############

# read in table of genes belonging to orthogroups, select significant OGs
# and parse into gene lists
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
GOtables_filenames <- list.files(pattern = "GOterm_mapping.tsv",  
                                 full.names = TRUE)[-c(4,5,10)]

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

# Parsing function (table to GeneSetCollection)
GeneSetCollection_constructor <- function(genome_GO_table = GOtable){
  #genome_GO_table=GO_tables[[3]]
  
  go_df_simple <- data.table(genome_GO_table %>%
                               dplyr::select(GO_ID, GeneID) %>%
                               filter(!is.na(GO_ID)))
  
  go_species <- genome_GO_table %>%
    dplyr::select(Genome) %>% 
    unique() %>% unlist() %>% unname()
  
  # Run once - takes ~ 1hr
  gmap <- buildGOmap(gomap = go_df_simple)
  saveRDS(gmap, file = paste("gmap.", go_species, ".rds", sep=""))
  
  full_universe <- genome_GO_table %>%
    dplyr::select(GeneID) %>%
    unique() %>% unlist() %>% unname()
  
  return(full_universe)
}

# Parse GO data table and make ready for clusterprofiler
GO_tables <- lapply(GOtables_filenames, GOtable_parsing)
full_gene_sets <- lapply(GO_tables, GeneSetCollection_constructor)
names(full_gene_sets) <- purrr::map(full_gene_sets, "Species")

# Enrichment analysis function
go.test.clusterpro <- function(full_set, gene_lists_tbl_test, namei) {
  # full_set=full_gene_sets[[3]]; gene_lists_tbl_test=gene_lists_tbl_contract; namei="contract"
  
  # Extract species ID
  go_species <- data.table(gene=full_set) %>%
    mutate(spp=tstrsplit(gene, "_")[[1]]) %>% 
    dplyr::select(spp) %>% unique() %>% unlist() %>% unname()
  
  print(go_species)
  
  # GO analysis
  gene.cand <- gene_lists_tbl_test %>% 
    dplyr::select(paste(go_species, ".protein",sep="")) %>% 
    unlist(use.names = F) %>% 
    strsplit(., split=",") %>%
    unlist() %>% 
    paste(go_species, "_", ., sep="") %>% 
    gsub(., pattern = " ", replacement = "") %>% 
    unlist() %>% 
    unique()
  
  # Universe - longest transcript proteins
  univ <- full_gene_sets[full_gene_sets %like% go_species][[1]] %>% unique() %>% unlist()
  
  spp.gmap = readRDS(list.files(path = "", 
                                pattern = paste("gmap.", go_species,".rds", sep=""), 
                                full.names = T))
  
  # Test enrichment
  go.test <- enricher(gene = gene.cand,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.01,
                      universe = univ, 
                      TERM2GENE = spp.gmap)
  
  # Find term description for each GO term
  go.test2 <- data.table(go.test@result) %>% filter(p.adjust<0.05)
  pp <- foreach(i=1:length(go.test2$ID), .combine = "rbind", .errorhandling = "remove") %do% {
    tmp.df = go.test2 %>%
      filter(ID == go.test2$ID[i])
    
    message(paste(i, go.test2$ID[i], sep = " | "))
    
    go.test2[ID==go.test2$ID[i]] %>% 
      mutate(term=go2term(ID)$Term) %>% 
      mutate(x = eval(parse(text=GeneRatio)),
             y = eval(parse(text=BgRatio))) %>% 
      mutate(odd = x/y)
  }
  
  plot.go <- {pp %>% 
      filter(p.adjust < 0.05) %>% 
      arrange(desc(p.adjust)) %>% 
      ggplot(aes(x = log2(odd),
                 y = reorder(term, log2(odd)),
                 size = Count,
                 color = log10(p.adjust))) +
      geom_point() +
      labs(x="log2(Enrichment)", 
           y="",
           color="FDR P-value",
           size="Candidate gene number") +
      theme_bw() +
      scale_color_viridis(option="magma") +
      #paletteer::scale_color_paletteer_c("gameofthrones::targaryen") +
      theme(title = element_text(face="bold", size=15),
            legend.text = element_text(face="bold", size=14),
            legend.title = element_text(face="bold", size=16),
            legend.background = element_blank(),
            axis.text.x = element_text(face="bold", size=15),
            axis.text.y = element_text(face="bold", size=15),
            axis.title.x = element_text(face="bold", size=18),
            axis.title.y = element_text(face="bold", size=18))}
  
  ggsave(plot.go, filename = paste(namei,".", go_species,".goterm.cafe.pdf", sep=""), width=20, height=12)
  
  # Write list for REVIGO 
  fwrite(x = pp %>% 
           select(ID, p.adjust, odd), 
         file = paste(namei,".", go_species,".goterm.cafe.revigo.txt", sep=""), 
         quote = F, row.names = F, col.names = F, sep = "\t")
  
  # Finish
  return(pp)
  
}

# Output enriched GO terms - contracted genes: go.out.con <- go.out
go.out.con <- lapply(c(full_gene_sets, gene_lists_tbl_contract, "contract"), 
                 FUN = go.test.clusterpro)

# Output enriched GO terms - expanded genes
go.out.exp <- lapply(c(full_gene_sets, gene_lists_tbl_expand, "expanded"), 
                 FUN = go.test.clusterpro)

# Output enriched GO terms for extractions & contractions
go.out.con <- lapply(con.list[-c(1,2,10)], go.test.clusterpro, set="contract")
go.out.exp <- lapply(ext.list[-c(1,2,10)], go.test.clusterpro, set="expand")

# Add GO term reduction
revigo_list_contract <- list.files(pattern = "contract.feb13.revigo", full.names = T)
revigo_list_expand <- list.files(pattern = "expand.feb13.revigo", full.names = T)

# Semantic reduction
library(rrvgo)

reduce.go <- function(go_list_enriched) {
  # gene_list_enriched=revigo_list_expand[[3]]
  
  # Extract species ID
  geneii=data.table(fread(gene_list_enriched))
  colnames(geneii) <- c("ID", "qvalue", "OR")
  go_species <- tstrsplit(gsub(x=gene_list_enriched, 
                               pattern="/scratch/csm6hg/data/",""),".", fixed=T)[[1]]
  print(go_species)
  
  # Conditional based on length
  if(nrow(geneii)<15) {
    print("List is very short, no redundancy analysis preformed.")
    reduced.data <- data.table(geneii, species=go_species) 
  } else {
    
    # Sim reduction
    simMatrix <- calculateSimMatrix(geneii$ID,
                                    orgdb="org.Dm.eg.db",
                                    ont="BP",
                                    method="Rel")
    
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
go.reduce.out.con <- lapply(go.out.con, reduce.go, mc.cores = threads)

# Finalize output
go.reduce.fin.con <- data.table(rbindlist(go.reduce.out.con, fill=T))

# Save output
saveRDS(go.reduce.out.con, file = "reduced.terms.contracting.rds")
saveRDS(go.reduce.out.con, file = "reduced.terms.expanding.rds")




### Rapidly evolving gene familes ###

# What are the most common GO terms in the rapidly evolving gene set?
rap <- list.files(path = "../data/",pattern = "sep9.revigo.txt", full.names = T)
rapi <- data.table(rbindlist(lapply(rap, fread), idcol=T))
colnames(rapi) <- c("spp", "GO", "p", "OR")
rapii <- data.table(rapi %>% 
          group_by(GO) %>% 
          count(GO))

# Rename species
rapi.o <- rapi %>%
  mutate(spp=case_when(spp==1 ~ "Carinata",
                        spp==2 ~ "Galeata",
                        spp==3 ~ "Magna",
                        spp==4 ~ "Pulex_euro",
                        spp==5 ~ "Pulex_nam",
                        spp==6 ~ "Pulicaria",
                        spp==7 ~ "Sinensis"))

# Write table data
#write.csv(rapi.o, file = "../data/rapidlyevolving.csv")