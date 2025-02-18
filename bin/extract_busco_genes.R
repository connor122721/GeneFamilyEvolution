#!/usr/bin/env Rscript
# Extract BUSCO genes across species - make fasta files for each species
# Connor Murray 7.30.2023, mod 1.20.2025
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(stringr)
library(foreach)
library(seqinr)
library(optparse)

# Parameters
option_list <- list(
  make_option(c("--max_missing_ingroup"), type="integer", default=0, help="Maximum number of ingroup species that can be missing a BUSCO gene"),
  make_option(c("--max_missing_outgroup"), type="integer", default=2, help="Maximum number of outgroup species that can be missing a BUSCO gene"),
  make_option(c("--metadata"), type="character", default="/project/berglandlab/connor/GeneFamilyEvolution/list.proteins", help="Metadata file"),
  make_option(c("--wd", "-w"), type="character", help="Path to output for nf-pipeline.")
)

opt <- parse_args(OptionParser(option_list=option_list))

max_missing_ingroup <- opt$max_missing_ingroup
max_missing_outgroup <- opt$max_missing_outgroup
metadata <- opt$metadata
wd <- opt$wd

print(max_missing_ingroup)
print(metadata)

# Testing
# max_missing_ingroup = 0; max_missing_outgroup = 0; metadata = "/project/berglandlab/connor/GeneFamilyEvolution/ncbi_genomes"; wd="/project/berglandlab/connor/GeneFamilyEvolution/output"

# Files
files.bus <- list.files(recursive = T, 
          path=wd, 
          pattern = "full_table", 
          full.names = T)

# Metadata
meta <- fread(metadata, header=F) %>% 
  mutate(species = V2, 
         group = V3,
         fasta = paste(wd, "/genomes/proteomes/", V2, ".protein.faa",sep="")) %>% 
  select(species, group, fasta)

# Read in busco files
l <- lapply(files.bus, fread, skip=2, fill=TRUE, sep="\t")
setattr(l, 'names', basename(str_remove(str_remove(files.bus, "full_table_"), ".tsv")))

# Bind all files
dt <- data.table(rbindlist(l, use.names = T, idcol = "species"))
if (ncol(dt) > 0) {
  colnames(dt)[c(2,7)] <- c("busco","url")
}

# Join with metadata
dt <- dt %>% left_join(meta, by="species")

# Filter for common single copy orthologs (SCO)
table(dt$Status, dt$species)

# Extract highest "Score" only when the Status is "Duplicated"
dt <- data.table(dt %>% group_by(species, busco) %>% filter(Status=="Complete"))

# Summarize the number of complete genes per busco and group
dt1 <- data.table(dt[Status == "Complete"] %>% group_by(busco, group) %>% summarize(n=n()))

# Filter for busco genes that are present in all ingroup species and allow up to two missing in the outgroups
ingroup_species <- unique(meta[group == "ingroup"]$species)
outgroup_species <- unique(meta[group == "outgroup"]$species)
sco <- dt1[busco %in% dt1[group == "ingroup" & n >= (length(ingroup_species) - max_missing_ingroup)]$busco & 
           busco %in% dt1[group == "outgroup" & n >= (length(outgroup_species) - max_missing_outgroup)]$busco]

dt2 <- dt[busco %in% sco$busco][Status == "Complete"] %>% arrange(busco)

# Extract and rename BUSCO complete genes for each species
foreach(i=1:length(unique(dt2$species))) %do% {

  spp <- unique(dt2$species)[i]
  print(spp)
  
  sequences <- dt2[species==spp]$Sequence
  busco_ids <- dt2[species==spp]$busco
  
  # Read the protein fasta file
  fasta_file <- meta[species == spp]$fasta
  fasta_data <- read.fasta(fasta_file, seqtype = "AA")
  
  # Filter the fasta data to include only the sequences in dt2
  filtered_fasta <- fasta_data[names(fasta_data) %in% sequences]
  
  # Ensure the sequences and BUSCO IDs are in the same order
  filtered_fasta <- filtered_fasta[order(match(names(filtered_fasta), sequences))]
  busco_ids <- busco_ids[order(match(sequences, names(filtered_fasta)))]
  
  # Rename the sequences to the BUSCO IDs
  names(filtered_fasta) <- busco_ids
  
  # Write the renamed sequences to a new fasta file
  write.fasta(sequences = filtered_fasta, 
    names = names(filtered_fasta), 
    file.out = paste(spp, ".rename.sco.faa", sep=""))
}

# Write unique BUSCO IDs to file
write.table(unique(dt2$busco), 
            file = "sco.txt", 
            quote = F, 
            row.names = F, 
            col.names = F)
