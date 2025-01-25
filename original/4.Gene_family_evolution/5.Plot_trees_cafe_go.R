# GOTerm plotting across expanding / contracting
# Connor Murray 2.26.2024
# ijob -c 15 --mem=50G -p standard -A berglandlab_standard
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(AnnotationDbi)
library(phytools)
library(ape)
library(tidyverse)
library(foreach)
library(ggtree)
library(viridis)

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

### Rapidly evolving gene familes ###

# What are the most common GO terms in the rapidly evolving gene set?
rap <- list.files(path = "../data",pattern = "_Revigo_BP_OnScreenTable.tsv", full.names = T)
rapi <- lapply(rap, fread, fill=T)
names(rapi)  <- rap
rapi <- data.table(rbindlist(rapi, fill=T, idcol=T, use.names =TRUE)) %>% 
  mutate(spp=gsub(x = tstrsplit(.id, "_")[[1]], pattern = "../data/", ""),
         mode=gsub(x = tstrsplit(.id, "_")[[2]], pattern = "../data/", ""))
colnames(rapi)[1:2] <- c("file", "GO")
rapii <- data.table(rapi %>% 
          group_by(GO) %>% 
          count(GO))

# Expanding
rapEx <- list.files(path = "../data",pattern = ".goterm.expand.feb13.revigo.txt", full.names = T)
rapiEx <- data.table(rbindlist(lapply(rapEx, fread, fill=T), 
                               fill=T, idcol=T, use.names =TRUE), mode="expand")
colnames(rapiEx) <- c("spp", "GO", "p", "OR", "mode")

# Contracting
rapCo <- list.files(path = "../data",pattern = ".goterm.contract.feb13.revigo.txt", full.names = T)
rapiCo <- data.table(rbindlist(lapply(rapCo, fread, fill=T), 
                               fill=T, idcol=T, use.names =TRUE), mode="contract")
colnames(rapiCo) <- c("spp", "GO", "p", "OR", "mode")

# Merge
dt1 <- data.table(rbind(rapiEx, rapiCo)) %>% 
  mutate(spp=case_when(spp==1 ~ "carinata",
                        spp==2 ~ "galeata",
                        spp==3 ~ "magna",
                        spp==4 ~ "pulex.euro",
                        spp==5 ~ "pulex.nam",
                        spp==6 ~ "pulicaria",
                        spp==7 ~ "sinensis"))

# Rename species
euro.rap <- dt1[spp=="pulex.euro"] %>%
  left_join(rapi, by=c("spp", "mode", "GO"))

rapi.o <- data.table(rbind(dt1 %>%
  left_join(rapi, by=c("spp", "mode", "GO")) %>% 
  filter(c(!is.na(file))),
  euro.rap) %>% 
  left_join(rapii, by="GO"))

# Find term description for each GO term
pp <- foreach(i=1:length(rapi.o$GO), .combine = "rbind", .errorhandling = "remove") %do% {
  print(i)
  data.table(rapi.o[i] %>% 
    mutate(term=clusterProfiler::go2term(GO)$Term))
}
 
# Write table data
# write.csv(rapi.o, file = "../data/expanding_contractingGO.csv")
# write.csv(pp, file = "../data/expanding_contractingGO_list_supp.csv")

# Extract top terms per species
ii <- data.table(rbind(pp[n>2], 
            pp[spp=="pulex.euro"]) %>%
    group_by(spp, mode) %>% 
    arrange(desc(-p)))
  

# Plot all go terms
plot.go1 <- {
  ii %>% 
    ggplot(., 
           aes(x = reorder(as.factor(term), 
                           log10(p)),
               y = -log10(p),
               fill = factor(spp),
               group = 1)) +
    geom_bar(position = position_dodge2(preserve = "single"), stat = "identity") +
    facet_wrap(~mode, nrow = 2) +
    labs(x = "", 
         y = "log10(FDR P-value)",
         fill = "Species") +
    theme_bw() +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=14),
          legend.title = element_text(face="bold", size=16),
          legend.background = element_blank(), 
          strip.text =  element_text(face="bold", size=15),
          axis.text.x = element_text(face="bold", size=15, angle=90),
          axis.text.y = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18))
}

#ggsave(plot = plot.go1, "/scratch/csm6hg/figs/ExpCon_GO_col.pdf", 
       dpi=300, height=12, width=10)

# Add overall groupings for Figure
ii <- data.table(ii %>% 
  mutate(grp = case_when(term %like% "sperm" | 
                         term %like% "male" | 
                         term %like% "germ" | 
                         term %like% "recom" |
                         term %like% "repro" ~ "Spermatogenesis and reproduction",
                         term %like% "iron" | 
                         term %like% "repair" | 
                         term %like% "oxido" ~ "Stress response", 
                         TRUE ~ "General metabolic processes"
                               )))
# Change order
ii$term <- factor(ii$term, levels = unique(ii$term[order(decreasing = F,ii$grp)]))

# Go term plot
plot.go2 <- {ii %>% 
    ggplot(., aes(x = as.factor(spp),
                  y = term,
                  size = log2(OR),
                  color = log10(p))) +
    geom_point() +
    facet_wrap(~mode, nrow=1) +
    labs(x="Species", 
         y="",
         color="log10(FDR P-value)",
         size="log2(Enrichment)") +
    theme_bw() +
    scale_color_viridis(option="magma") +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=14),
          legend.title = element_text(face="bold", size=16),
          legend.background = element_blank(),
          strip.text = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=15, angle=-45),
          axis.text.y = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18))}

ggsave(plot = plot.go2, "/scratch/csm6hg/figs/ExpCon_GO_new_sperm.pdf", dpi=300)

# VennDiagram of GO terms
library(eulerr)

# Create a list of sets
iii <- data.table(pp[mode=="expand"])
setsExp <- list(
  Euro.Dpulex = iii[spp=="pulex.euro"]$GO,
  NAm.Dpulex = iii[spp=="pulex.nam"]$GO,
  Dpulicaria = iii[spp=="pulicaria"]$GO,
  Dgaleata = iii[spp=="galeata"]$GO,
  Dmagna = iii[spp=="magna"]$GO,
  Dsinensis = iii[spp=="sinensis"]$GO,
  Dcarinata = iii[spp=="carinata"]$GO
)

iiii <- data.table(pp[mode=="contract"])
setsCon <- list(
  Euro.Dpulex = iiii[spp=="pulex.euro"]$GO,
  NAm.Dpulex = iiii[spp=="pulex.nam"]$GO,
  Dpulicaria = iiii[spp=="pulicaria"]$GO,
  Dgaleata = iiii[spp=="galeata"]$GO,
  Dmagna = iiii[spp=="magna"]$GO,
  Dsinensis = iiii[spp=="sinensis"]$GO,
  Dcarinata = iiii[spp=="carinata"]$GO
)

# Create the Venn diagram
vennExp <- euler(setsExp)
vennCon <- euler(setsCon)

# Plot the Venn diagram
pdf("/scratch/csm6hg/figs/Exp_GO_Venn_All.pdf", width=8, height=6)
plot(vennExp, quantities = TRUE, lty = 1:3, labels = list(font = 4))
dev.off()

pdf("/scratch/csm6hg/figs/Con_GO_Venn_All.pdf", width=8, height=6)
plot(vennCon, quantities = TRUE, lty = 1:3, labels = list(font = 4))
dev.off()