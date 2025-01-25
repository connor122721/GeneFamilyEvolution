# absrel output and odds ratio of expanding and contracting genes undergoing selection
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)

# Working directory 
setwd("/scratch/csm6hg/phd/codeml/")

# Hog -> OG metadata
og <- data.table(fread("/scratch/csm6hg/phd/cafe_sep7/N0.sep7.orthofinder.tsv"))

# Hogs annotation
hog.ann <- data.table(fread("/scratch/csm6hg/phd/cafe_sep7/hogs.sep7.function.tsv") %>% 
                        left_join(og, by=c("Gene"="HOG")))

# Parsed aBSREL data
dt1 <- data.table(readRDS("/scratch/csm6hg/phd/data/parsed.meta.absrel.ogs.rds"))

# Read in contrasting terms
con <- data.table(readRDS("/scratch/csm6hg/phd/data/contracted.genes.sep7.rds"))

# Read in expanding terms
exp <- data.table(readRDS("/scratch/csm6hg/phd/data/expanded.genes.sep7.rds"))

# Average omega across exp vs cont.
dt2 <- dt1[omega1<10][!species=="node"] %>% 
  mutate(evo=case_when(Gene %in% exp[Gene_Tree=="n0"]$tree ~ "Expanding",
                       Gene %in% con[Gene_Tree=="n0"]$tree ~ "Contracting",
                       TRUE ~ "Non-Fluctuating Genes")) 

# How many species present in tree tested
hognum <- data.table(fread("/scratch/csm6hg/phd/cafe_sep7/hog_sep7_gene_counts.tsv"))

# Restrict to genes in Cafe result 
dt.num <- data.table(dt2[Gene %in% unique(hognum$HOG), {
  blank_row_count <- sum(rowSums(is.na(.SD) | .SD == "") > 0)
  .(blank_row_count)}, 
  by = .(OG, file, evo), 
  .SDcols=colnames(dt2 %>% select(contains("protein")))])
OGsKeep <- unique(dt.num[blank_row_count<9]$OG)

# Restrict to genes without too many overrepresented expansions; Genes=5
num=5
HOGsNumKeep <- unique(hognum[!c(carinata>num|pulex.nam>num|pulex.euro>num|
                    pulicaria>num|sinensis>num|magna>num|galeata>num)]$HOG)

# Which genes are significantly evolving
dt3 <- data.table(dt2[OG %in% OGsKeep][Gene %in% HOGsNumKeep][Gene_Tree=="n0"] %>% 
  group_by(evo, OG) %>%
  summarize(sig=case_when(adj_pval<0.05 & 
         c(omega1>1 & omega1<10 | omega2>1 & omega2<10 | omega3>1 & omega3<10) ~"Y",
                          TRUE ~ "N")) %>% 
  group_by(evo, OG) %>% 
  distinct() %>% 
  group_by(OG) %>% 
  mutate(OG_count=n()) %>% 
  group_by(evo,OG) %>% 
  summarize(sig=case_when(OG_count==2 ~ "Y",
                          TRUE ~ sig)))

dt3 <- dt3[!duplicated(dt3)]

# Odds ratio of expansion and contraction
or.exp <- dt3 %>%
  filter(evo != "Contracting") %>%
  select(evo, sig) %>%
  table() %>%
  `[`(, c("Y", "N")) %>%  # Reorder the columns
  fisher.test()

or.con <- dt3 %>%
  filter(evo != "Expanding") %>%
  select(evo, sig) %>%
  table() %>%
  `[`(, c("Y", "N")) %>%  # Reorder the columns
  fisher.test()

# Odds ratios
or.con
or.exp

# Stats
t.test(as.numeric(dt2[evo=="Contracting"]$omega1), 
       as.numeric(dt2[evo=="Expanding"]$omega1))
