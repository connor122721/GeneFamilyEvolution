
# Investigation of gene family evolution among _Daphnia_
This repository contains a comprehensive analysis of gene family dynamics within several species of _Daphnia_, freshwater microcrustaceans that serve as valuable model systems in evolutionary genetics. This study examines how gene family expansion and contraction contribute to ecological adaptation and evolutionary innovation across _Daphnia_ species. We apply comparative genomics approaches to explore the forces driving gene family evolution, focusing on key functional categories like spermatogenesis and stress response genes, and evaluate selection pressures shaping these gene families. Our findings provide insights into the interplay between gene family turnover and adaptive processes, demonstrating both general patterns and species-specific dynamics in gene family evolution within _Daphnia_.

## BUSCO
BUSCO (Benchmarking Universal Single-Copy Orthologs) is used to assess the completeness of the _Daphnia_ genome assemblies and annotations by evaluating conserved single-copy orthologs. 

## Orthologs
Identifying orthologous genes across species is essential for understanding evolutionary relationships.
_OrthoFinder_ was used to detect orthologs within and across _Daphnia_ species.

## Phylogenomics
Phylogenomic analysis was used to infer evolutionary relationships and dynamics within _Daphnia_ species. 
We used _MCMCtree_ on BUSCO genes.

## Gene family evolution
This section explores gene family expansion and contraction across _Daphnia_ species, focusing on the genetic basis of ecological adaptations. 
We used _Cafe5_ and _ClusterProfiler_.
Expansion of genes related to spermatogenesis and stress responses across taxa, suggesting adaptation to environmental pressures.

## Selection 
Selection analysis investigates evolutionary pressures on specific gene families, particularly those undergoing expansion, such as spermatogenesis and stress response genes.
General trends of positive selection in expanding gene families and species-specific patterns of purifying selection.
Use of codon-based models (e.g., _PAML_, _HyPhy_) to detect adaptive versus purifying selection.
