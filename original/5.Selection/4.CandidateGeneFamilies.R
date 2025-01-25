# Find candidate genes
# Connor Murray 3.25.2024
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(ggtree)

# Working directory
setwd("/scratch/csm6hg/cafe_sep7/")

# Common theme element
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

# Hog -> OG metadata
og <- data.table(fread("N0.sep7.orthofinder.tsv"))

# Hogs annotation
hog.ann <- data.table(fread("hogs.sep7.function.tsv") %>% 
                        left_join(og, by=c("Gene"="HOG")))

# ABSREL output
dt1 <- data.table(readRDS("/scratch/csm6hg/data/parsed.meta.absrel.ogs.rds"))

# Read in data
contract <- data.table(readRDS("/scratch/csm6hg/data/contracted.genes.sep7.rds"))
expand <- data.table(readRDS("/scratch/csm6hg/data/expanded.genes.sep7.rds"))

# Spermatogenesis terms
expand[OG %in% hog.ann[Annotation %like% "sperm"]$OG]
contract[OG %in% hog.ann[Annotation %like% "sperm"]$OG]
dt1[OG %in% hog.ann[Annotation %like% "sperm"]$OG][!branch %like% "Node"][!branch %like% "NODE"]

# Spermatogenesis terms
expand[OG %in% hog.ann[Annotation %like% "ferritin"|Annotation %like% "iron"|Annotation %like% "heme"]$OG]
contract[OG %in% hog.ann[Annotation %like% "ferritin"|Annotation %like% "iron"|Annotation %like% "heme"]$OG]
iron <- dt1[OG %in% hog.ann[Annotation %like% "ferritin"|Annotation %like% "iron"|Annotation %like% "heme"]$OG][!branch %like% "Node"][!branch %like% "NODE"]

iron[adj_pval<0.05]
