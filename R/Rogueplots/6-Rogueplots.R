################################################################################
# Name: 6-Rogueplots.R
# Authors: Lucas Buffan and Frédéric Lozes-Méléro
# Emails: lucas.l.buffan@gmail.com, frederic.lozes-melero@umontpellier.fr 
# Goal: Rogue plots to see wich fossils must be removed 
#       (Klopfstein & Spasojevic, 2019)
################################################################################

library(ape)

## Guet the function written by Klopfstein & Spasojevic (2019) -----------------
source("./R/Rogueplots/function_Klopfstein_2018.R")

## ConsensusTree file in nexus format ------------------------------------------
consensus.tree <- read.nexus("Path/to/consensus/tree")

## Different Run trees files in nexus format -----------------------------------
trees1 <- read.nexus("Path/to/run1/trees")
trees2 <- read.nexus("Path/to/run2/trees")

## Set the burnin fraction use in the analyses --------------------------------- 
burnin_frac <- 0.25
trees1 <- trees1[(length(trees1) * burnin_frac + 1):length(trees1)]
trees2 <- trees2[(length(trees2) * burnin_frac + 1):length(trees2)]
reftrees <- c(trees1, trees2)

## Set the fossils list for the rogue plots ------------------------------------ 
fossils_list = consensus.tree[["tip.label"]][1:58]

## Rogueplot command------------------------------------------------------------ 
for (x in 1:length(fossils_list)) {
  taxon <- fossils_list[x]
  plot <- paste0("./Results/Phylogeny/Rogueplots/Rogue_plot_", taxon, ".pdf")
  table <- paste0("./Results/Phylogeny/Rogueplots/Rogue_placement_", taxon, ".txt")
  create.rogue.plot(consensus.tree,
                    reftrees,
                    rogues = taxon,
                    outgroup = fossils_list[56], 
                    type = "Spectral",
                    col = NULL,
                    min.prob = 0.01,
                    outfile.table = table,
                    outfile.plot = plot)
}

## References------------------------------------------------------------------- 
# Klopfstein, S. & Spasojevic, T. (2019) Illustrating phylogenetic placement of fossils using RoguePlots: An example from ichneumonid parasitoid wasps (Hymenoptera, Ichneumonidae) and an extensive morphological matrix. PLOS ONE, 14, e0212942.


