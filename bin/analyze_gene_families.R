# Find ingroup vs outgroup interesting gene lists
# Connor Murray 1.23.2024; modified 2.1.2025
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(viridis)
library(argparse)

# Working directory
setwd("/project/berglandlab/connor/GeneFamilyEvolution/test")

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

# Read in expanded vs contracted genes
dt <- data.table(readRDS("expanded.genes.rds"))

# Read in metadata
meta <- data.table(fread("../list.proteins", header = F))
colnames(meta) <- c("proteome", "group")
meta[, spp := gsub("\\..*", "", basename(proteome))]

# Filter by ingroup status and tips of gene tree (most recent)
dt.filt <- dt %>% 
  left_join(meta, by = "spp") %>%
  filter(Gene_Tree == "n0") %>%
  group_by(tree) %>%
  filter(all(group == "ingroup"))

# Find common
dt.filt$Annotation
