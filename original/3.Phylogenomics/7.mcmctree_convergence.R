# Analyze MCMCtree output convergence from 2 runs/plot data
# code adapted from https://github.com/mariodosreis/divtime/blob/master/R/analysis.R#L37
# Connor Murray 
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(coda)
library(tidyverse)
library(MCMCtreeR)
library(patchwork)

# Working directory
setwd("C:/Users/csm6hg/Desktop/GeneFamilyEvolution/3.Phylogenomics/mcmctree_data/")

# Phylogenetic tree
phy <- readMCMCtree("FigTree_run1.tre", from.file = TRUE)
#ape::write.tree(phy$apePhy, file="mcmctree_busco_daphnia.nwk")

# Plot w/ time
MCMC.tree.plot(phy, 
               analysis.type = "MCMCtree", 
               cex.tips = 0.9,
               time.correction = 100, 
               relative.height = 0.08, 
               plot.type = "cladogram", 
               lwd.bar = 2, 
               label.offset = 4,
               scale.res = c("Eon", "Period"),
               node.method = "node.length", 
               col.age = "#008b0080",
               cex.labels = 1)

# Posterior times - reading data0 run and actual run
prior <- fread("mcmc_run_data0.txt")
mcmc1 <- fread("mcmc_run1.txt")
mcmc2 <- fread("mcmc_run2.txt")

# to check for convergence of the MCMC runs, we calculate the posterior
# means of times for each run, and plot them against each other
t.mean1 <- apply(mcmc1[,2:7], 2, mean)
t.mean2 <- apply(mcmc2[,2:7], 2, mean)

# Plot divergence times between runs
tt <- cbind(data.table(t.mean1), 
            data.table(t.mean2))

# Convergence 
con <- {
  
  tt %>% 
  ggplot(aes(x=t.mean1*100,
             y=t.mean2*100)) +
  geom_abline(slope = 1, intercept = 0, size=1) +
  geom_point(size=3, color="blue") +
  coord_fixed(ratio=1) +
  theme_bw() +
  labs(x="Mean Divergence Run 1 (mya)",
       y="Mean Divergence Run 2 (mya)") +
  theme(strip.text = element_text(face="bold", size=18),
        title = element_text(face="bold", size=20),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))
}

# Trace plots
gens <- seq(0, 5000000, by=10000)

# Merge runs
mcmc.merge <- data.table(rbind(data.table(mcmc1, run="Run_1"),
                               data.table(mcmc2, run="Run_2")))

tplot <-{
  mcmc.merge[Gen %in% gens] %>%
    select(Gen,run,colnames(mcmc1)[colnames(mcmc1) %like% "t_n"]) %>% 
    pivot_longer(cols = colnames(mcmc1)[colnames(mcmc1) %like% "t_n"]) %>% 
    ggplot(., aes(x=Gen/1e6, 
                  y=value*100, 
                  color=run)) +
      geom_line(size=1.2, alpha=0.7) +
      scale_color_manual(values=c("Run_1"="blue", "Run_2"="orange")) +
      facet_wrap(~name) +
      theme_bw() +
      labs(x="MCMC Generation (Million)",
           y="Divergence time (mya)",
           color="") +
      theme(strip.text = element_text(face="bold", size=18),
            title = element_text(face="bold", size=20),
            legend.text = element_text(face="bold", size=18),
            legend.position = "bottom",
            axis.text.x = element_text(face="bold", size=18),
            axis.text.y = element_text(face="bold", size=18),
            axis.title.x = element_text(face="bold", size=20),
            axis.title.y = element_text(face="bold", size=20),
            axis.title = element_text(face="bold", size=20))
}

# Save composite plot
comp <- (tplot | con) + 
  plot_layout(widths = c(2.5, 1)) +
  plot_annotation(tag_levels = 'A', tag_suffix = ".")

ggsave(plot=comp, filename = "mcmc.con1.pdf", dpi = 300, width = 16, height = 6)

# we can calculate the effective sample sizes (ESS) of the parameters
# (you need to have the coda package installed for this to work)
mean.mcmc <- apply(mcmc1[,-1], 2, mean)
ess.mcmc <- apply(mcmc1[,-1], 2, coda::effectiveSize)
var.mcmc <- apply(mcmc1[,-1], 2, var)
se.mcmc <- sqrt(var.mcmc / ess.mcmc)
cbind(mean.mcmc, ess.mcmc, var.mcmc, se.mcmc)

mean.mcmc <- apply(mcmc2[,-1], 2, mean)
ess.mcmc <- apply(mcmc2[,-1], 2, coda::effectiveSize)
var.mcmc <- apply(mcmc2[,-1], 2, var)
se.mcmc <- sqrt(var.mcmc / ess.mcmc)
cbind(mean.mcmc, ess.mcmc, var.mcmc, se.mcmc)
