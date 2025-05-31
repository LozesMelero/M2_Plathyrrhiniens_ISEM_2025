################################################################################
# Name: 9-BioGeoBEARS_DEC+J.R
# Authors: Lucas Buffan and Frédéric Lozes-Méléro
# Emails: lucas.l.buffan@gmail.com, frederic.lozes-melero@orange.fr 
# Goal: NWM's historical biogeography analysis using BioGeoBEARS (Matzke, 2013)
################################################################################

rm(list = ls())

library(ape)
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
# seems to sometimes fail on simple problems (2-3 parameters)
library(rexpokit)
library(cladoRcpp)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)

#sp <-  "all"
#sp <- "Haplorrhini"
sp <- "Anthropoidea"

All = "All_Primates/" ; Hap = "Haplorrhini/" ; Anth = "Anthropoidea/"

if (sp == "all") {
  phy <- paste("./Data/Historical_biogeography/Trees/All_extant/Cons_all_tree.newick")
  temp <- read.tree(phy)
  table_geo = paste("./Data/Historical_biogeography/relaxed_model/",All,"Area_biogeography_platy.txt",sep="")
  time_period = "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_time_period.txt"
  connectivity_matrix <- "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_area_matrix.txt"
  dispersal_matrix <- "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_dispersal_matrix.txt"
  
} else if (sp == "Haplorrhini") {
  phy <- paste("./Data/Historical_biogeography/Trees/All_extant/Haplo_cons_tree.newick")
  temp <- read.tree(phy)
  table_geo = paste("./Data/Historical_biogeography/relaxed_model/",Hap,"Area_biogeography_platy.txt",sep="")
  time_period = "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_time_period.txt"
  connectivity_matrix <- "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_area_matrix.txt"
  dispersal_matrix <- "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_dispersal_matrix.txt"
  
} else if (sp == "Anthropoidea") {
  phy <- paste("./Data/Historical_biogeography/Trees/All_extant/Anth_cons_tree.newick")
  temp <- read.tree(phy)
  table_geo = paste("./Data/Historical_biogeography/relaxed_model/",Anth,"Area_biogeography_platy.txt",sep="")
  time_period = "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_time_period_Anth.txt"
  connectivity_matrix <- "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_area_matrix_Anth.txt"
  dispersal_matrix <- "./Data/Historical_biogeography/relaxed_model/BioGeoBEARS_dispersal_matrix_Anth.txt"
} 


#### DEC+J analysis -------------------------------------------------------------
DECJr_platyrrhini = define_BioGeoBEARS_run()
DECJr_platyrrhini$trfn = phy 
DECJr_platyrrhini$geogfn = table_geo 

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=table_geo) ## Get the maximum range size observed
max(rowSums(dfnums_to_numeric(tipranges@df))) 
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))
DECJr_platyrrhini$max_range_size = max_range_size 

areas = getareas_from_tipranges_object(tipranges) ## Get area (letter)
states_area_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=T)
DECJr_platyrrhini$states_list = states_area_list 

DECJr_platyrrhini$timesfn = time_period 
DECJr_platyrrhini$timeperiods = unlist((unname(c(read.table(time_period)))))
#DECJr_platyrrhini$areas_adjacency_fn = connectivity_matrix 
DECJr_platyrrhini$dispersal_multipliers_fn = dispersal_matrix 
#DECJr_platyrrhini$areas_allowed_fn = "./Data/Historical_biogeography/BioGeoBEARS_allowedarea_matrix.txt"

DECJr_platyrrhini$min_branchlength = 0.000001    
DECJr_platyrrhini$include_null_range = TRUE
DECJr_platyrrhini$force_sparse = FALSE
DECJr_platyrrhini$on_NaN_error = -1e50
DECJr_platyrrhini$speedup = FALSE        
DECJr_platyrrhini$use_optimx = TRUE    
DECJr_platyrrhini$num_cores_to_use = 20

DECJr_platyrrhini = section_the_tree(inputs=DECJr_platyrrhini, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)

DECJr_platyrrhini$return_condlikes_table = TRUE
DECJr_platyrrhini$calc_TTL_loglike_from_condlikes_table = TRUE
DECJr_platyrrhini$calc_ancprobs = TRUE  

DECJr_platyrrhini = readfiles_BioGeoBEARS_run(DECJr_platyrrhini) 

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
if (sp == "all") {
  load("./Results/Historical_biogeography/All_Primates/DEC/Relaxed_model/Platyrrhini_DECR_allprimate.Rdata")
  resDECR = resR
} else if (sp == "Haplorrhini") {
  load("./Results/Historical_biogeography/Haplorrhini/DEC/Relaxed_model/Platyrrhini_DECR_Haplo.Rdata")  
  resDECR = resR
} else if (sp == "Anthropoidea") {
  load("./Results/Historical_biogeography/Anthropoidea/DEC/Relaxed_model/Platyrrhini_DECR_Anth.Rdata")  
  resDECR = resR
} 
dstart = resDECR$outputs@params_table["d","est"]
estart = resDECR$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
DECJr_platyrrhini$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DECJr_platyrrhini$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DECJr_platyrrhini$BioGeoBEARS_model_object@params_table["e","init"] = estart
DECJr_platyrrhini$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
DECJr_platyrrhini$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DECJr_platyrrhini$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DECJr_platyrrhini$BioGeoBEARS_model_object@params_table["j","est"] = jstart

DECJr_platyrrhini = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=DECJr_platyrrhini)
check_BioGeoBEARS_run(DECJr_platyrrhini)

if (sp == "all") {
  resfnjr = paste("./Results/Historical_biogeography/",All,"DECJ/Relaxed_model/Platyrrhini_DECJr_allprimate.Rdata",sep="")
} else if (sp == "Haplorrhini") {
  resfnjr = paste("./Results/Historical_biogeography/",Hap,"DECJ/Relaxed_model/Platyrrhini_DECJr_Haplo.Rdata",sep="")
} else if (sp == "Anthropoidea") {
  resfnjr = paste("./Results/Historical_biogeography/",Anth,"DECJ/Relaxed_model/Platyrrhini_DECJr_Anth.Rdata",sep="")
} 

runslow = F
if (runslow){
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  resjr = bears_optim_run(DECJr_platyrrhini)
  resjr    
  save(resjr, file=resfnjr)
  resDECjr = resjr
} else {
  load(resfnjr)
  resDECjr = resjr
}

#######################################################
#### Plot ancestral states - DECJ -----------------------------------------------
#######################################################
if (sp == "all") {
  pdffn = paste("./Results/Historical_biogeography/",All,"DECJ/Relaxed_model/Platyrrhini_DECJr_allprimate.pdf",sep="")
} else if (sp == "Haplorrhini") {
  pdffn = paste("./Results/Historical_biogeography/",Hap,"DECJ/Relaxed_model/Platyrrhini_DECJr_Haplo.pdf",sep="")
} else if (sp == "Anthropoidea") {
  pdffn = paste("./Results/Historical_biogeography/",Anth,"DECJ/Relaxed_model/Platyrrhini_DECJr_Anth.pdf",sep="")
} 
pdf(pdffn, height=30, width=7)
analysis_titletxt ="BioGeoBEARS DEC+J on platyrrhini w/ relaxed dispersal matrix"
tr = read.tree(phy)
results_object = resDECjr
scriptdir ="./R/Historical_Biogeography/a_scripts/"
# States
res2 = plot_BioGeoBEARS_results(results_object, addl_params=list("j"), plotwhat="text", label.offset=0.2, tipcex=0.35, statecex=0.25, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.2, tipcex=0.35, statecex=0.10, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

####################################
# Probabilities of states/ranges at each node
####################################
# To get the probabilities of each state/range at each node:
resDECjr$ML_marginal_prob_each_state_at_branch_top_AT_node

tr = read.tree(phy)
trtable = prt(tr, printflag=FALSE)
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))
areas = getareas_from_tipranges_object(tipranges)

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based)){    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}
ranges_list

# Make the node numbers the row names
# Make the range_list the column names
range_probabilities = as.data.frame(resjr$ML_marginal_prob_each_state_at_branch_top_AT_node)
row.names(range_probabilities) = trtable$node
names(range_probabilities) = ranges_list
range_probabilities = round(range_probabilities, 2)

if (sp == "all") {
  write.table(range_probabilities, file=paste("./Results/Historical_biogeography/",All,"DECJ/Relaxed_model/range_probabilities_primates.txt",sep=""), quote=FALSE, sep="\t")
} else if (sp == "Haplorrhini") {
  write.table(range_probabilities, file=paste("./Results/Historical_biogeography/",Hap,"DECJ/Relaxed_model/range_probabilities_Haplo.txt",sep=""), quote=FALSE, sep="\t")
} else if (sp == "Anthropoidea") {
  write.table(range_probabilities, file=paste("./Results/Historical_biogeography/",Anth,"DECJ/Relaxed_model/range_probabilities_Anth.txt",sep=""), quote=FALSE, sep="\t")
} 

if (sp == "all") {
  nodes_of_interest <- c(243, 248, 414, 249, 459, 456, 450)
} else if (sp == "Haplorrhini") {
  nodes_of_interest <- c(227, 232, 398, 233, 443, 440, 432)
} else if (sp == "Anthropoidea") {
  nodes_of_interest <- c(227, 393, 228, 438, 435, 427)
} 

subset_range_probs <- range_probabilities[as.character(nodes_of_interest), ]

if (length(nodes_of_interest) == 7) {
  new_names <- c("Haplorrhini",
                 "Catarrhini/Platyrrhini",
                 "Catarrhini",
                 "Platyrrhini",
                 "Ashaninkacebus/spp",
                 "Perupithecus/spp",
                 "Ucayalipithecus/spp")
} else {
  new_names <- c("Catarrhini/Platyrrhini",
                 "Catarrhini",
                 "Platyrrhini",
                 "Ashaninkacebus/spp",
                 "Perupithecus/spp",
                 "Ucayalipithecus/spp")
} 

rownames(subset_range_probs) <- new_names

if (sp == "all") {
  write.table(subset_range_probs, file=paste("./Results/Historical_biogeography/",All,"DECJ/Relaxed_model/subset_range_probabilities_primates.txt",sep=""), quote=FALSE, sep="\t")
} else if (sp == "Haplorrhini") {
  write.table(subset_range_probs, file=paste("./Results/Historical_biogeography/",Hap,"DECJ/Relaxed_model/subset_range_probabilities_Haplo.txt",sep=""), quote=FALSE, sep="\t")
} else if (sp == "Anthropoidea") {
  write.table(subset_range_probs, file=paste("./Results/Historical_biogeography/",Anth,"DECJ/Relaxed_model/subset_range_probabilities_Anth.txt",sep=""), quote=FALSE, sep="\t")
} 

####################################
# Calculate per-area probabilities for each node
####################################
relprobs_matrix = resDECjr$ML_marginal_prob_each_state_at_branch_top_AT_node
probs_each_area = round(infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based),2)

colnames(probs_each_area) <- areas
rownames(probs_each_area) <- trtable$node
subset_probs_each_area <- probs_each_area[as.character(nodes_of_interest), ]
rownames(subset_probs_each_area) <- new_names

if (sp == "all") {
  write.table(subset_probs_each_area, file=paste("./Results/Historical_biogeography/",All,"DECJ/Relaxed_model/per-area_probabilities_primates.txt",sep=""),
              sep ="\t", row.names = T, col.names = T, quote = FALSE)
} else if (sp == "Haplorrhini") {
  write.table(subset_probs_each_area, file=paste("./Results/Historical_biogeography/",Hap,"DECJ/Relaxed_model/per-area_probabilities_Haplo.txt",sep=""),
              sep ="\t", row.names = T, col.names = T, quote = FALSE)
} else if (sp == "Anthropoidea") {
  write.table (subset_probs_each_area, file=paste("./Results/Historical_biogeography/",Anth,"DECJ/Relaxed_model/per-area_probabilities_Anth.txt",sep=""),
               sep ="\t", row.names = T, col.names = T, quote = FALSE)
} 

## References ------------------------------------------------------------------
# Matzke, N.J. (2013) BioGeoBEARS: BioGeography with Bayesian (and likelihood) evolutionary analysis in R Scripts. R package, version 0.2, 1, 2013.