# Extract longest ortholog for European D. pulex (exception genome)
# Connor Murray 7.30.2023, mod 1.20.2025
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(stringr)
library(seqinr)

# Read in GFF
gff <- data.table(fread("/project/berglandlab/connor/genomes/pulex_euro/Daphnia.aed.0.6.gff") %>% 
                  mutate(len=V5-V4,
                         gene=str_remove(tstrsplit(V9, ";")[[1]], "ID="),
                         parent_gene=str_remove(tstrsplit(V9, ";")[[2]], "Parent=")))

# Restrict to protein-coding genes
gff <- gff[V3=="mRNA"]

# Select largest transcript
gff.dt <- data.table(gff %>% group_by(parent_gene) %>% top_n(1, len))

# Still duplicated
dup <- unique(gff.dt[duplicated(gff.dt$parent_gene)]$parent_gene)

# Sliced genes
dup.genes <- data.table(gff.dt[parent_gene %in% dup] %>% group_by(parent_gene) %>% sample_n(1))$gene

# Fin primary transcripts
fin <- unique(gff.dt[!parent_gene %in% dup]$gene,
              dup.genes)

# Output primary transcripts
write.table(fin, "euro_primary_transcripts.txt", quote = F, row.names = F, col.names = F)

# Extract primary transcripts from pulexeuro.protein.faa
pulex_fasta <- read.fasta("/project/berglandlab/connor/genomes/proteins_species/pulexeuro.protein.faa", seqtype = "AA")
primary_transcripts <- pulex_fasta[names(pulex_fasta) %in% fin]

# Write primary transcripts to a new fasta file
write.fasta(sequences = primary_transcripts, 
            names = names(primary_transcripts), 
            file.out = "pulexeuro.protein.faa")