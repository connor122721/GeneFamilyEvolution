# Extract BUSCO genes across species
# Connor Murray 7.30.2023

# Libraries
library(data.table)
library(tidyverse)
library(stringr)
library(foreach)
library(phylotools)

# Working directory
setwd("/scratch/csm6hg/genomes/proteins_species/primary_transcripts")

# Files
files.bus <- list.files(recursive = T, pattern = "full_table.tsv")[-c(4,5,10)]

# Read in busco files
l <- lapply(files.bus, fread, skip=2, fill=TRUE, sep="\t")
setattr(l, 'names', tstrsplit(files.bus, "/")[[1]])

# Bind all files
dt <- data.table(rbindlist(l, use.names = T, idcol = T))
colnames(dt)[c(1,2,7)] <- c("species","busco","url")

# Filter for common single copy orthologs (SCO)
table(dt$Status, dt$species)
dt1 <- data.table(dt[Status=="Complete"] %>% group_by(busco) %>% summarize(n=n()))
sco <- dt1[n==length(unique(dt$species))]
dt2 <- dt[busco %in% sco$busco][Status=="Complete"] %>% arrange(busco)

# write SCO for each species
foreach(i=1:length(unique(dt2$species))) %do% {

  spp <- unique(dt2$species)[i]
  print(spp)
  
  write.table(dt2[species==spp]$Sequence, 
              file = paste(spp, ".sco.txt", sep=""), 
              quote = F, 
              row.names = F, 
              col.names = F)
  
  ref_table <- data.frame(og=dt2[species==spp]$Sequence,
                          new=dt2[species==spp]$busco)
  
  rename.fasta(infile = paste(spp, ".sco.faa", sep=""), 
               ref_table, 
               outfile = paste(spp, ".rename.sco.faa", sep=""))
  
}

write.table(unique(dt2$busco), 
            file = "sco.txt", 
            quote = F, 
            row.names = F, 
            col.names = F)
