# Extract longest ortholog for European D. pulex (exception genome)
# Connor Murray 7.30.2023, mod 1.20.2025
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(stringr)
library(seqinr)
library(readxl)

# Read in GFF
gff <- data.table(fread("/project/berglandlab/connor/GeneFamilyEvolution/output/genomes/gtf/pulexeuro.gff") %>% 
                  mutate(len=V5-V4,
                         gene=str_remove(tstrsplit(V9, ";")[[1]], "ID="),
                         parent_gene=str_remove(tstrsplit(V9, ";")[[2]], "Parent=")))

# Read in gene annotations
panth <- data.table(read_excel("/project/berglandlab/connor/GeneFamilyEvolution/output/genomes/gtf/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# Primary sequences
prime <- c(fread("pulexeuro.sequences", header=F))

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
pulex_fasta <- read.fasta("/project/berglandlab/connor/GeneFamilyEvolution/output/genomes/proteomes/pulexeuro.protein.faa", seqtype = "AA")
primary_transcripts <- pulex_fasta[names(pulex_fasta) %in% fin]

# Write primary transcripts to a new fasta file
write.fasta(sequences = primary_transcripts, 
            names = names(primary_transcripts), 
            file.out = "pulexeuro_orf.protein.faa")

# Rename fasta headers
library(phylotools)
rename.fasta(infile = "pulexeuro_orf.protein.faa", 
             ref_table = panth2,
             outfile = "pulexeuro_rename.protein.faa")
