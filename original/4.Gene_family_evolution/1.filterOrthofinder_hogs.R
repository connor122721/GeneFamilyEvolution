# Filter HOGs for CAFE5
# Connor Murray 9.7.2023
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

library(data.table)
library(tidyverse)

setwd("/scratch/csm6hg/")

ortho <- data.table(fread("genomes/primary/OrthoFinder/Results_Sep07/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"))

# Filter large gene counts in one species and exclude genes present in only 1 species
hog <- ortho[, OG := NULL]
hog[, `Gene Tree Parent Clade` := NULL]
hog <- melt(hog, id.vars='HOG', variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

# Exclude HOGs with lots of genes in a one or more species. 
# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max < 100]$HOG
hog <- hog[HOG %in% keep]

# Exclude HOGs present in only 1 species
keep <- hog[, .N, HOG][N > 1]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'genomes/primary/hog_sep7_gene_counts.tsv', sep='\t')
system("sed -i 's/.protein//g' genomes/primary/hog_sep7_gene_counts.tsv")

