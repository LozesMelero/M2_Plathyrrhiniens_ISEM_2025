################################################################################
# Name: DNA-sequences-taxa.R
# Author: Frédéric Lozes-Méléro
# Email: frederic.lozes-melero@umontpellier.fr 
# Goal: Extract a list of platyrrhini species with and without DNA sequences 
# output from Genophy
################################################################################

## Raw list from Genophy Output ------------------------------------------------
list <- read.delim("./Data/Taxonomy/cOX3.fasta", header = F)

## List of taxa from Genophy raw Output ----------------------------------------
species_genophy<- list[grep(">", list[,1]), ]
species_genophy <- gsub("^>", "", species_genophy) # remove ">"

## List of undefined species/taxa from Genophy Output --------------------------
toMatch <- c("_sp.", "_cf.", "_X_", "_x_")
undefined_species <- species_genophy[grep(paste(toMatch,collapse="|"), species_genophy)]

species_genophy_kept <- species_genophy[species_genophy %in% undefined_species == F]

## Concatenate taxonomy (handbook) ---------------------------------------------
platy_list <- read.delim("./Data/Taxonomy/Full_species_list_Lynx2024.txt", header = T)
# Genus species subspecies (e.g. Alouatta_palliata_palliata)
platy_genus_species_subspecies <- ifelse(
  platy_list$Subspecies == "", 
  paste(platy_list$Genus, platy_list$Species, sep = "_"),
  paste(platy_list$Genus, platy_list$Species, platy_list$Subspecies, sep = "_"))
  # Genus species (for backtracing, e.g. Alouatta_geoffroyi)
platy_genus_species <- paste(platy_list$Genus, platy_list$Species, sep = "_")

## Get absent species ----------------------------------------------------------
  # Round 1: full species not in taxo list
ABS_species_round1 <- setdiff(species_genophy_kept, platy_genus_species_subspecies)
  # Round 2: subspecies does not exist
ABS_species_round2 <- setdiff(ABS_species_round1, platy_genus_species)
  # Round 3: subspecies now turned to species
ABS_genus_subspecies <- sapply(X = ABS_species_round2,
                                      FUN = function(x){
                                        spl = strsplit(x, split = "_")[[1]]
                                        # if name of length 3, i.e. subspecies
                                        if(length(spl) == 3){
                                          return(paste0(spl[1], "_", spl[3]))
                                          }
                                        else{
                                          return(NA)
                                        }})
ABS_genus_subspecies <- unname(unlist(ABS_genus_subspecies))
missing <- which(ABS_genus_subspecies %in% platy_genus_species_subspecies)
ABS_species_round3 <- ABS_species_round2[-missing]
saveRDS(ABS_species_round3, "./Data/missing_in_genophy.RDS")
