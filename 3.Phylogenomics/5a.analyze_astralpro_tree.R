# Extract BUSCO genes across species
# Connor Murray 7.30.2023
# ; R

library(ape)
library(ggtree)
library(phytools)
library(data.table)

# Working directory
setwd("C:/Users/Conno/Desktop/gitmaster/chapter2_code/Gene_family/3.Phylogenomics/")

# Species-tree from Astral-pro
tre <- read.tree("sco_species_astralpro")

# BUSCO gene trees from IQtree2
trees <- read.tree("busco_sco_genes.tre")[-c(266)]

# Check integrity
test <- rbindlist(lapply(X = trees, function(x) {
  data.table(num.taxa=length(unique(x$tip.label)),
             class=class(x),
             nnode=x$Nnode)}))

# Set root to carinata
tre.root <- root.phylo(tre, outgroup = "carinata")
trees.root <- root.multiPhylo(trees, outgroup = "carinata")

# ggtree(force(tree = tre.root, method = "extend")) + geom_nodelab() + geom_tiplab()
# ggtree(trees.root[4]) + geom_nodelab() + geom_tiplab()

# Make tree ultrametric for plotting
trees.root2 <- lapply(trees.root, FUN = function(x) {force.ultrametric(x)})
class(trees.root2) <- "multiPhylo"
trees.root2 <- root.multiPhylo(trees.root2, outgroup = "carinata", resolve.root=T)

is.ultrametric.multiPhylo(phy = trees.root2)
is.binary.multiPhylo(phy = trees.root2)
is.rooted.multiPhylo(phy = trees.root2)

class(trees.root2)
# ggtree(trees.root2[[4]]) + geom_nodelab() + geom_tiplab()

# Plot ultrametric trees - densitrees
density1 <- trees.root2 %>% 
  ggdensitree(branch.length='none',
              alpha=0.1,
              align.tips = TRUE,
              colour='steelblue') + 
  geom_tiplab(size=5)

# Save
ggsave(plot = density1, 
       filename = "densitree.daphnia.ultrametric.sco.new.pdf", 
       width=10, 
       height=10)

# Output
# write.tree(trees.root2, file = "busco_sco_genes_rooted_ultrametric.tre")

# Plot - species tree
p <- tre.root %>% 
  ggtree(branch.length='none') + 
  geom_nodelab() +
  geom_tiplab(size=3)

# Save
ggsave(plot = p, 
       filename = "astralpro.sco.pdf", 
       width=10, 
       height=10)
