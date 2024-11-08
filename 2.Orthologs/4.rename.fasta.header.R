# Add functional annotations to D.pulex euro protein fasta
# 8.13.2023
# ijob -c 1 --mem=50G -p largemem -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(readxl)

# Working directory
setwd("/scratch/csm6hg/genomes/primary")

# Read in gene annotations
panth <- data.table(read_excel("/project/berglandlab/daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Primary sequences
# system("cat pulex.euro.protein.faa.fai | cut -f1 > pulex.euro.sequences")
prime <- c(fread("pulex.euro.sequences", header=F))

# restrict to primary sequences
panth1 <- panth[qseqid %in% prime$V1]

# Some transcripts are not annotated: 605
unk <- data.table(qseqid=prime$V1[which(!prime$V1 %in% panth1$qseqid)],
                  pro=paste(prime$V1[which(!prime$V1 %in% panth1$qseqid)], 
                            "Uncharacterized protein [Daphnia pulex euro]", 
                            sep=" "))

# New header for each sequence: [old header, new header]
panth2 <- data.table(rbind(panth1 %>% 
  group_by(qseqid) %>% 
  summarize(pro = paste(qseqid,
                        protein_name, 
                        "[Daphnia pulex euro]", 
                        sep=" ")),
  unk))

#write.table(panth2, "pulex.euro/pulex.euro.protein.rename.txt", sep="\t", quote = F, row.names = F, col.names = F)

# Rename fasta headers
library(phylotools)
rename.fasta(infile = "pulex.euro.protein.sub.faa", 
             ref_table = panth2,
             outfile = "pulex.euro.protein.rename.faa")
