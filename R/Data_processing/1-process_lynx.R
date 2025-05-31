################################################################################
# Name: 1-process_lynx.R
# Authors: Lucas Buffan and Frédéric Lozes-Méléro
# Emails: lucas.l.buffan@gmail.com, frederic.lozes-melero@umontpellier.fr 
# Goal: List of platyrrhine species from the book "Neotropical Primates" 
#       Published by Lynx Nature Books in association with Re:wild in 2024
################################################################################

library(tidyverse)
source("./R/useful/helper_functions.R")

## Raw list from Lynx Nature Books (2024) Neotropical Primates -----------------
raw <- read.delim("./Data/Taxonomy/Lynx_checklist_raw.txt")

# IUCN Red List Categories (acronyms)
l_iucn <- c("EX", "EW", "RE", "CR", "EN", "NE", "VU", "NT", "LC", "DD" , "NA")

# Genera list (24) of New World Monkey (following the book)   
list_genera <- c( "Aotus", "Saguinus", "Leontocebus", "Oedipomidas", "Tamarinus", "Cebuella", 
                  "Callibella","Callimico", "Mico", "Callithrix", "Leontopithecus", "Cebus", 
                  "Sapajus", "Saimiri", "Alouatta", "Ateles", "Lagothrix", "Brachyteles", 
                  "Chiropotes", "Cacajao", "Pithecia", "Callicebus", "Plecturocebus", "Cheracebus")

## Extracting and create species list ------------------------------------------
species_list <- data.frame(species = character())

 # Code for extracting species name by using IUCN acronyms and Genus name
for (i in 1:nrow(raw)) {
  txt <- raw[i,1]
  spl <- strsplit(txt, split = " ")[[1]]
  index_iucn <- which(spl %in% l_iucn)
  
  if (length(index_iucn) > 0) {
    if(spl[index_iucn - 2] %in% list_genera){ 
      name <- paste(spl[index_iucn - 2], spl[index_iucn - 1])  
      species_list <- rbind(species_list, data.frame(species = name))
    }
    else if(spl[index_iucn - 3] %in% list_genera){  
      name <- paste(spl[index_iucn - 3],spl[index_iucn - 2], spl[index_iucn - 1])  
      species_list <- rbind(species_list, data.frame(species = name))
    }
 }
}

## Three-columns table ---------------------------------------------------------
  # Function to extract the desired feature from a species name
get_taxo <- function(x, what = c("genus", "species", "subspecies")){
  if(what == "genus"){
    i <- 1
  }
  else if(what == "species"){
    i <- 2
  }
  else if(what == "subspecies"){
    i <- 3
  }
  extracted <- strsplit(x, split = " ")[[1]][i]
  return(extracted)
}
  # Modify `species_list` so it now has 3 columns (gen, sp, subsp)
spl_new <- species_list %>% 
  mutate(Genus = sapply(X = species, FUN = get_taxo, what = "genus"),
         Species = sapply(X = species, FUN = get_taxo, what = "species"),
         Subspecies = sapply(X = species, FUN = get_taxo, what = "subspecies")) %>%
  select(!(species))
  # Save
write.tbl.std(spl_new, "./Data/Taxonomy/Full_species_list_Lynx2024.txt")

## Some references -------------------------------------------------------------
# Lynx Nature Books (2024) Neotropical Primates. Lynx Illustrated Checklists. Lynx Nature Books, Barcelona.
