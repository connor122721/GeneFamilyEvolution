# GFF - get largest transcript per gene
# Connor Murray 7.28.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(stringr)

# Working directory
setwd("/scratch/csm6hg/genomes/")

# Read in GFF
gff <- data.table(fread("pulex_euro/Daphnia.aed.0.6.gff") %>% 
                  mutate(len=V5-V4,
                         gene=str_remove(tstrsplit(V9, ";")[[1]], "ID="),
                         parent_gene=str_remove(tstrsplit(V9, ";")[[2]], "Parent=")) %>% 
                  mutate(iso=str_remove(tstrsplit(gene, "-")[[2]], "R")))

# Restrict to protein-coding genes
gff <- gff[V3=="mRNA"]

# Select largest transcript - pick by alphabet for ties
gff.dt <- data.table(gff %>% group_by(parent_gene) %>% top_n(1, len) %>% arrange(iso) %>% top_n(1))

# Fin primary transcripts
fin <- unique(gff.dt$gene)

# output primary transcripts
write.table("pulex_euro/primary_transcripts.txt", x=fin, quote = F, row.names = F, col.names = F)
