# Extract BUSCO genes across species
# Connor Murray 7.30.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

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
              file = paste("/scratch/csm6hg/genomes/busco/", spp, ".sco.txt", sep=""), 
              quote = F, 
              row.names = F, 
              col.names = F)
  
  ref_table <- data.frame(og=dt2[species==spp]$Sequence,
                          new=dt2[species==spp]$busco)
  
  rename.fasta(infile = paste("/scratch/csm6hg/genomes/sco_aug28/", spp, ".sco.faa", sep=""), 
               ref_table, 
               outfile = paste("/scratch/csm6hg/genomes/sco_aug28/rename/", spp, ".rename.sco.faa", sep=""))
  
}

write.table(unique(dt2$busco), 
            file = "/scratch/csm6hg/genomes/sco/sco.aug28.txt", 
            quote = F, 
            row.names = F, 
            col.names = F)


# Past this point - only specific to Daphnia pule euro

# Read in GFF
gff <- data.table(fread("pulex_euro/Daphnia.aed.0.6.gff") %>% 
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

# output primary transcripts
write.table("pulex_euro/primary_transcripts.txt", x=fin, quote = F, row.names = F, col.names = F)
