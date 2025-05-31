################################################################################
# Name: 5-Molecular_backbone.R
# Authors: Lucas Buffan and Frédéric Lozes-Méléro
# Emails: lucas.l.buffan@gmail.com, frederic.lozes-melero@umontpellier.fr 
# Goal: Establish partial constraints for total evidence dating analysis, 
#       to facilitate convergence of markov chains
################################################################################

library(ape)
library(paleotree)

## Tree file in nexus format ---------------------------------------------------
tree <-read.nexus("./Results/Phylogeny/Molecular_phylogeny/Platy_molphylogeny.nex")

## Create partial constraints file ---------------------------------------------
createMrBayesConstraints(tree, partial = T, file = "./Results/Phylogeny/Molecular_phylogeny/Backbone/TED_Partial_constraint.txt")


# Bapst, D.W. (2012) paleotree: an R package for paleontological and phylogenetic analyses of evolution. Methods in Ecology and Evolution, 3, 803–807.
