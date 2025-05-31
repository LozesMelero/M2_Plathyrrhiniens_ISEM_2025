################################################################################
# Name: 7-Biogeography-data-processing.R
# Authors: Lucas Buffan and Frédéric Lozes-Méléro
# Emails: lucas.l.buffan@gmail.com, frederic.lozes-melero@orange.fr 
# Goal: Establish adjacency and dispersal matrix for historical biogeography 
#       analysis using BioGeoBEARS (Matzke, 2013)
################################################################################

library("tidygraph")
library("igraph")
library("ggraph")
library("tidyverse")
library("ape")

## Species repartition file in tsv format --------------------------------------
table_biogeo <- read.table("./Data/Historical_biogeography/Area_biogeography_platy.tsv", sep ="\t", header = TRUE)

## Area connectivity trough defined time period in tsv format ------------------
table <- read.table("./Data/Historical_biogeography/Connectivity_through_time.tsv", sep ="\t", header = TRUE)

## ConsensusTree file in nexus format ------------------------------------------
phy <- read.nexus("./Data/Historical_biogeography/Trees/All_extant/Platy_TED_FBD_datation_3clock_10rogueless_40M.con.tre")

## Remove unwanted species
#sp <-  "all"
#sp <- "Haplorrhini"
sp <- "Anthropoidea"

if (sp == "all") {
  print("All primates included in the phylogeny")
  
} else if (sp == "Haplorrhini") {
  rm_tip <- c(phy$tip.label[1:11],phy$tip.label[238:240])
  phy <- drop.tip(phy,rm_tip,rooted = is.rooted(phy))
  write.tree(phy, file = "./Data/Historical_biogeography/Trees/All_extant/Haplo_cons_tree.newick")
  print("Only Haplorrhini are included in the phylogeny")
  
} else if (sp == "Anthropoidea") {
  rm_tip <- c(phy$tip.label[1:11],phy$tip.label[234:240])
  phy <- drop.tip(phy,rm_tip,rooted = is.rooted(phy))
  write.tree(phy, file = "./Data/Historical_biogeography/Trees/All_extant/Anth_cons_tree.newick")
  print("Only Anthropoidea are included in the phylogeny")
} 

## Keeping only species present in consensusTree file (.tre) -------------------
## (e.g. useful for rogueless and/or reduce species matrix analysis)------------
  table_biogeo <- table_biogeo[gsub(" ", "_", table_biogeo$Taxa) %in% phy$tip.label,]

## Creating final table for BioGeoBEARS analysis (phylip format) ---------------
## and exporting .txt file (needed for BioGeoBEARS script) ---------------------
table_DEC_biogeo <- rbind(
  c(length(phy$tip.label), ncol(table_biogeo)-1),  # Number of area 
  cbind(table_biogeo$Taxa,as.matrix(apply(table_biogeo[, 2:ncol(table_biogeo)], 1, paste, collapse = "")))
)

All = "All_primates/"
Hap = "Haplorrhini/"
Anth = "Anthropoidea/"

if (sp == "all") {
  write.table(table_DEC_biogeo, paste("./Data/Historical_biogeography/constrained_model/",All,"Area_biogeography_platy.txt",sep=""), 
              sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(table_DEC_biogeo, paste("./Data/Historical_biogeography/relaxed_model/",All,"Area_biogeography_platy.txt",sep=""), 
              sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
} else if (sp == "Haplorrhini") {
  write.table(table_DEC_biogeo, paste("./Data/Historical_biogeography/constrained_model/",Hap,"Area_biogeography_platy.txt",sep=""), 
              sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(table_DEC_biogeo, paste("./Data/Historical_biogeography/relaxed_model/",Hap,"Area_biogeography_platy.txt",sep=""), 
              sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
} else if (sp == "Anthropoidea") {
  write.table(table_DEC_biogeo, paste("./Data/Historical_biogeography/constrained_model/",Anth,"Area_biogeography_platy.txt",sep=""), 
              sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(table_DEC_biogeo, paste("./Data/Historical_biogeography/relaxed_model/",Anth,"Area_biogeography_platy.txt",sep=""), 
              sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
} 

## This function create both adjacency and dispersal matrix using the ----------
## connectivity trough period/time file ----------------------------------------
## Adjacency matrix <- line 47 to 60 -------------------------------------------
## First, a squared matix, for each time period is created fill with "0" -------
## Then, this script replace all "0" with a "1" if areas are connected according 
## to the connectivity trough time table ---------------------------------------
## Dispersal matrix <- 67 to 105 -----------------------------------------------
## Lets consider 4 areas: A,B,C & D. Lets also consider that A & B are connected
## [AB], then B & C [BC] and then C & D [CD] -----------------------------------
## First, a squared matix, for each time period is created fill with "0.0000001"
## This represent the dispersal rate between 2 separate areas ------------------
## Then, this script replace all concerned "0.0000001" with either : -----------
## 1 if this is the same area [A->A], [B->B],-----------------------------------
## 0.5 if 2 areas are connected [A->B], [B->C], [C->D] -------------------------
## 0.25 if 2 areas are not connected but an area connects them [A->C], [A->B->C] 
## 0.125 if 2 areas are not connected but 2 area connects them [A->D], [A->B->C->D] 
## Then all the matrix are exported to .txt file, use in BioGeoBEARS script ----
make_table_BioGeoBEARS <- function(time_periods_table, prefix){
  temp_length <- (ncol(time_periods_table) - 2)
  table_adjacency <- c()
  table_dispersal <- c()
  temp_mat_adj <- matrix(0, temp_length, temp_length)
  for(i in 1:nrow(time_periods_table)){
    temp_mat_adj <- matrix(0, temp_length, temp_length)
    temp_period_table <- time_periods_table[i, c(3:ncol(time_periods_table))]
    from <- c()
    to <- c()
    
    for(j in 1:temp_length){
      if(temp_period_table[j] != "0"){
        temp_mat_adj[j,eval(parse(text = temp_period_table[j]))] <- 1
        from <- c(from, rep(j,length(eval(parse(text = temp_period_table[j])))))
        to <- c(to, eval(parse(text = temp_period_table[j])))
      }
    }
    
    nodes <- tibble(id = 1:temp_length)
    edges <- tibble(from = from, to = to)
    
    temp_data_from <- edges[edges[,1] != edges[,2],]
    temp_data_to_1 <- edges[edges[,1] != edges[,2],]
    temp_data_to_2 <- edges[edges[,1] != edges[,2],]
    
    colnames(temp_data_to_1) <- c("to", "to_2")
    colnames(temp_data_to_2) <- c("to_2", "to_3")
    
    temp_data_03 <- merge(temp_data_from, temp_data_to_1, by = "to")
    temp_data_03 <- temp_data_03[temp_data_03[,2] != temp_data_03[,3],]
    
    temp_data_04 <- merge(temp_data_03, temp_data_to_1, by = "to_2")
    temp_data_04 <- temp_data_04[temp_data_04[,3] != temp_data_04[,4],]
    
    temp_mat_dispersal <- matrix(0.0000001, ncol(table_biogeo)-1, ncol(table_biogeo)-1)
    temp_data_direct <- as.matrix(edges)
    
    for(k in 1:nrow(temp_data_direct)){
      if(temp_data_direct[k,1] == temp_data_direct[k,2]){
        temp_mat_dispersal[temp_data_direct[k,1], temp_data_direct[k,2]] <- 1 
      }
      if(temp_data_direct[k,1] != temp_data_direct[k,2]){
        temp_mat_dispersal[temp_data_direct[k,1], temp_data_direct[k,2]] <- 0.5
      }
    }
    
    for(k in 1:nrow(temp_data_03)){
      if(temp_mat_dispersal[temp_data_03[k,2], temp_data_03[k,3]] == 0.0000001){
        temp_mat_dispersal[temp_data_03[k,2], temp_data_03[k,3]] <- 0.25
      }
    }
    
    for(k in 1:nrow(temp_data_04)){
      if(temp_mat_dispersal[temp_data_04[k,3], temp_data_04[k,4]] == 0.0000001){
        temp_mat_dispersal[temp_data_04[k,3], temp_data_04[k,4]] <- 0.125
      }
    }
    
    table_adjacency <- rbind(table_adjacency, LETTERS[1:ncol(table_biogeo)-1], temp_mat_adj, matrix(data=" ", ncol=ncol(table_biogeo)-1, nrow=1))
    table_dispersal <- rbind(table_dispersal, LETTERS[1:ncol(table_biogeo)-1], temp_mat_dispersal, matrix(data=" ", ncol=ncol(table_biogeo)-1, nrow=1))
  }
  table_adjacency <- rbind(table_adjacency, cbind("END", matrix(data = " ", ncol=ncol(table_adjacency)-1, nrow=1)))
  table_dispersal <- rbind(table_dispersal, cbind("END", matrix(data = " ", ncol=ncol(table_adjacency)-1, nrow=1)))
  write.table(time_periods_table[,2], paste(prefix, "_time_period.txt", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(table_adjacency, paste(prefix, "_area_matrix.txt", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(table_dispersal, paste(prefix, "_dispersal_matrix.txt", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
data <- make_table_BioGeoBEARS(table, "./Data/Historical_biogeography/constrained_model/BioGeoBEARS")

## Relaxed model with 1, 0.75, 0.5, 0.25, 0.1 values ---------------------------
make_table_BioGeoBEARS_relaxed <- function(time_periods_table, prefix){
  temp_length <- (ncol(time_periods_table) - 2)
  table_adjacency <- c()
  table_dispersal <- c()
  temp_mat_adj <- matrix(0, temp_length, temp_length)
  for(i in 1:nrow(time_periods_table)){
    temp_mat_adj <- matrix(0, temp_length, temp_length)
    temp_period_table <- time_periods_table[i, c(3:ncol(time_periods_table))]
    from <- c()
    to <- c()
    
    for(j in 1:temp_length){
      if(temp_period_table[j] != "0"){
        temp_mat_adj[j,eval(parse(text = temp_period_table[j]))] <- 1
        from <- c(from, rep(j,length(eval(parse(text = temp_period_table[j])))))
        to <- c(to, eval(parse(text = temp_period_table[j])))
      }
    }
    
    nodes <- tibble(id = 1:temp_length)
    edges <- tibble(from = from, to = to)
    
    temp_data_from <- edges[edges[,1] != edges[,2],]
    temp_data_to_1 <- edges[edges[,1] != edges[,2],]
    temp_data_to_2 <- edges[edges[,1] != edges[,2],]
    
    colnames(temp_data_to_1) <- c("to", "to_2")
    colnames(temp_data_to_2) <- c("to_2", "to_3")
    
    temp_data_03 <- merge(temp_data_from, temp_data_to_1, by = "to")
    temp_data_03 <- temp_data_03[temp_data_03[,2] != temp_data_03[,3],]
    
    temp_data_04 <- merge(temp_data_03, temp_data_to_1, by = "to_2")
    temp_data_04 <- temp_data_04[temp_data_04[,3] != temp_data_04[,4],]
    
    temp_mat_dispersal <- matrix(0.1, ncol(table_biogeo)-1, ncol(table_biogeo)-1)
    temp_data_direct <- as.matrix(edges)
    
    for(k in 1:nrow(temp_data_direct)){
      if(temp_data_direct[k,1] == temp_data_direct[k,2]){
        temp_mat_dispersal[temp_data_direct[k,1], temp_data_direct[k,2]] <- 1 
      }
      if(temp_data_direct[k,1] != temp_data_direct[k,2]){
        temp_mat_dispersal[temp_data_direct[k,1], temp_data_direct[k,2]] <- 0.75
      }
    }
    
    for(k in 1:nrow(temp_data_03)){
      if(temp_mat_dispersal[temp_data_03[k,2], temp_data_03[k,3]] == 0.1){
        temp_mat_dispersal[temp_data_03[k,2], temp_data_03[k,3]] <- 0.5
      }
    }
    
    for(k in 1:nrow(temp_data_04)){
      if(temp_mat_dispersal[temp_data_04[k,3], temp_data_04[k,4]] == 0.1){
        temp_mat_dispersal[temp_data_04[k,3], temp_data_04[k,4]] <- 0.25
      }
    }
    
    table_adjacency <- rbind(table_adjacency, LETTERS[1:ncol(table_biogeo)-1], temp_mat_adj, matrix(data=" ", ncol=ncol(table_biogeo)-1, nrow=1))
    table_dispersal <- rbind(table_dispersal, LETTERS[1:ncol(table_biogeo)-1], temp_mat_dispersal, matrix(data=" ", ncol=ncol(table_biogeo)-1, nrow=1))
  }
  table_adjacency <- rbind(table_adjacency, cbind("END", matrix(data = " ", ncol=ncol(table_adjacency)-1, nrow=1)))
  table_dispersal <- rbind(table_dispersal, cbind("END", matrix(data = " ", ncol=ncol(table_adjacency)-1, nrow=1)))
  write.table(time_periods_table[,2], paste(prefix, "_time_period.txt", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(table_adjacency, paste(prefix, "_area_matrix.txt", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(table_dispersal, paste(prefix, "_dispersal_matrix.txt", sep = ""), sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
data <- make_table_BioGeoBEARS_relaxed(table, "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS")

## References ------------------------------------------------------------------
# Matzke, N.J. (2013) BioGeoBEARS: BioGeography with Bayesian (and likelihood) evolutionary analysis in R Scripts. R package, version 0.2, 1, 2013.

