################################################################################
# Name: 8-BioGeoBEARS_DEC_DECJ.R
# Authors: Lucas Buffan and Frédéric Lozes-Méléro
# Emails: lucas.l.buffan@gmail.com, frederic.lozes-melero@orange.fr 
# Goal: NWM's historical biogeography analysis using BioGeoBEARS (Matzke, 2013)
################################################################################

library(ape)
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
# seems to sometimes fail on simple problems (2-3 parameters)
library(rexpokit)
library(cladoRcpp)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)

## ConsensusTree file in nexus format ------------------------------------------
phy <- paste("./Data/Historical_biogeography/Trees/All_extant/Platy_TED_FBD_datation_3clock_10rogueless_40M.con.tre.newick")
temp <- read.tree(phy)

## Import file generated in data-processing R Script ---------------------------
## Species repartition file in txt format --------------------------------------
table_geo = "./Data/Historical_biogeography/Area_biogeography_platy.txt"

## Defined time period ---------------------------------------------------------
time_period = "./Data/Historical_biogeography/BioGeoBEARS_time_period.txt"
#time_period = "./Data/Historical_biogeography/period.txt"

## Area connectivity trough defined time period (matrix for each time period) --
connectivity_matrix <- "./Data/Historical_biogeography/BioGeoBEARS_area_matrix.txt"

## Dispersal matrix trough defined time period (matrix for each time period) ---
dispersal_matrix <- "./Data/Historical_biogeography/BioGeoBEARS_dispersal_matrix.txt"

#### DEC analysis --------------------------------------------------------------
DEC_platyrrhini = define_BioGeoBEARS_run() ## Intitialize the model
DEC_platyrrhini$trfn = phy ## Location of the phylogeny file
DEC_platyrrhini$geogfn = table_geo ## Location of the geography text file

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=table_geo) ## Get the maximum range size observed
max(rowSums(dfnums_to_numeric(tipranges@df))) # here -> 3
max_range_size = 3
DEC_platyrrhini$max_range_size = max_range_size # Set the maximum number of areas any species may occupy

areas = getareas_from_tipranges_object(tipranges) ## Get area (letter)
states_area_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=T)
DEC_platyrrhini$states_list = states_area_list ## All possible state combination

# Set up a time-stratified analysis:
DEC_platyrrhini$timesfn = time_period ## Location of the time period file
DEC_platyrrhini$timeperiods = unlist((unname(c(read.table(time_period)))))
#DEC_platyrrhini$areas_adjacency_fn = connectivity_matrix ## adjacency matrix
DEC_platyrrhini$dispersal_multipliers_fn = dispersal_matrix ## Dispersal matrix
#DEC_platyrrhini$areas_allowed_fn = "./Data/Historical_biogeography/BioGeoBEARS_allowedarea_matrix.txt"

DEC_platyrrhini$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DEC_platyrrhini$include_null_range = TRUE

DEC_platyrrhini$force_sparse = FALSE
DEC_platyrrhini$on_NaN_error = -1e50
DEC_platyrrhini$speedup = FALSE        
DEC_platyrrhini$use_optimx = TRUE    
DEC_platyrrhini$num_cores_to_use = 8

#Divide the tree up by timeperiods/strata
DEC_platyrrhini = section_the_tree(inputs=DEC_platyrrhini, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)

DEC_platyrrhini$return_condlikes_table = TRUE
DEC_platyrrhini$calc_TTL_loglike_from_condlikes_table = TRUE
DEC_platyrrhini$calc_ancprobs = TRUE    # get ancestral states from optim run

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
DEC_platyrrhini = readfiles_BioGeoBEARS_run(DEC_platyrrhini) 

# checking inputs
DEC_platyrrhini = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=DEC_platyrrhini)
check_BioGeoBEARS_run(DEC_platyrrhini)

runslow = T
resfn = "./Results/Historical_biogeography/DEC/Normal_model/Platyrrhini_DEC_v1.Rdata"
if (runslow){
  res = bears_optim_run(DEC_platyrrhini)
  res    
  save(res, file=resfn)
  resDEC = res
} else {
  load(resfn)
  resDEC = res
}

#######################################################
#### Plot ancestral states - DEC -----------------------------------------------
#######################################################
pdffn = "./Results/Historical_biogeography/DEC/Normal_model/Platyrrhini_DEC_disp.pdf"
pdf(pdffn, height=15, width=7)
analysis_titletxt ="BioGeoBEARS DEC on platyrrhini w/ dispersal matrix"
tr = read.tree(phy)
results_object = resDEC
scriptdir ="./R/Historical_Biogeography/a_scripts/"
# States
res2 = plot_BioGeoBEARS_results(results_object, addl_params=list("j"), plotwhat="text", label.offset=0.2, tipcex=0.35, statecex=0.25, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.2, tipcex=0.35, statecex=0.10, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#-------------------------------------------------------------------------------

#### DEC+J analysis -------------------------------------------------------------

DECJ_platyrrhini = define_BioGeoBEARS_run()
DECJ_platyrrhini$trfn = phy 
DECJ_platyrrhini$geogfn = table_geo 
DECJ_platyrrhini$max_range_size = max_range_size 
DECJ_platyrrhini$states_list = states_area_list 
DECJ_platyrrhini$timesfn = time_period 
DECJ_platyrrhini$timeperiods = unlist((unname(c(read.table(time_period)))))
#DECJ_platyrrhini$areas_adjacency_fn = connectivity_matrix 
DECJ_platyrrhini$dispersal_multipliers_fn = dispersal_matrix 
#DECJ_platyrrhini$areas_allowed_fn = "./Data/Historical_biogeography/BioGeoBEARS_allowedarea_matrix.txt"
DECJ_platyrrhini$min_branchlength = 0.000001    
DECJ_platyrrhini$include_null_range = TRUE
DECJ_platyrrhini$force_sparse = FALSE
DECJ_platyrrhini$on_NaN_error = -1e50
DECJ_platyrrhini$speedup = FALSE        
DECJ_platyrrhini$use_optimx = TRUE    
DECJ_platyrrhini$num_cores_to_use = 8
DECJ_platyrrhini = section_the_tree(inputs=DECJ_platyrrhini, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
DECJ_platyrrhini$return_condlikes_table = TRUE
DECJ_platyrrhini$calc_TTL_loglike_from_condlikes_table = TRUE
DECJ_platyrrhini$calc_ancprobs = TRUE   
DECJ_platyrrhini = readfiles_BioGeoBEARS_run(DECJ_platyrrhini) 

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
DECJ_platyrrhini$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DECJ_platyrrhini$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DECJ_platyrrhini$BioGeoBEARS_model_object@params_table["e","init"] = estart
DECJ_platyrrhini$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
DECJ_platyrrhini$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DECJ_platyrrhini$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DECJ_platyrrhini$BioGeoBEARS_model_object@params_table["j","est"] = jstart

DECJ_platyrrhini = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=DECJ_platyrrhini)
check_BioGeoBEARS_run(DECJ_platyrrhini)

resfn = "./Results/Historical_biogeography/DECJ/Normal/Platyrrhini_DECJ_v1.Rdata"
runslow = T
if (runslow){
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  res = bears_optim_run(DECJ_platyrrhini)
  res    
  save(res, file=resfn)
  resDECj = res
} else {
  load(resfn)
  resDECj = res
}

#######################################################
#### Plot ancestral states - DECJ -----------------------------------------------
#######################################################
pdffn = "./Results/Historical_biogeography/DECJ/Normal/Platyrrhini_DECJ_disp.pdf"
pdf(pdffn, height=15, width=7)
analysis_titletxt ="BioGeoBEARS DEC+J on platyrrhini w/ dispersal matrix"
tr = read.tree(phy)
results_object = resDECj
scriptdir ="./R/Historical_Biogeography/a_scripts/"
# States
res2 = plot_BioGeoBEARS_results(results_object, addl_params=list("j"), plotwhat="text", label.offset=0.2, tipcex=0.35, statecex=0.25, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.2, tipcex=0.35, statecex=0.10, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#-------------------------------------------------------------------------------

#### DEC null model analysis ---------------------------------------------------
DEC_platyrrhini_null = define_BioGeoBEARS_run()
DEC_platyrrhini_null$trfn = phy
DEC_platyrrhini_null$geogfn = table_geo
DEC_platyrrhini_null$max_range_size = max_range_size
DEC_platyrrhini_null$min_branchlength = 0.000001   
DEC_platyrrhini_null$include_null_range = TRUE   
DEC_platyrrhini_null$on_NaN_error = -1e50  
DEC_platyrrhini_null$speedup = TRUE         
DEC_platyrrhini_null$use_optimx = TRUE   
DEC_platyrrhini_null$num_cores_to_use = 8
DEC_platyrrhini_null$force_sparse = FALSE  
DEC_platyrrhini_null$return_condlikes_table = TRUE
DEC_platyrrhini_null$calc_TTL_loglike_from_condlikes_table = TRUE
DEC_platyrrhini_null$calc_ancprobs = TRUE    

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=DEC_platyrrhini_null)
check_BioGeoBEARS_run(DEC_platyrrhini_null)

runslow = F
resfn = "./Results/Historical_biogeography/DEC/Null_model/Platyrrhini_DEC_v1.Rdata"
if (runslow){
  res = bears_optim_run(DEC_platyrrhini_null)
  res    
  save(res, file=resfn)
  resDECnull = res
} else {
  load(resfn)
  resDECnull = res
}

#######################################################
#### Plot ancestral states - DEC null Model ------------------------------------
#######################################################
pdffn = "./Results/Historical_biogeography/DEC/Null_model/Platyrrhini_DEC_null_model.pdf"
pdf(pdffn, height=15, width=7)
analysis_titletxt ="BioGeoBEARS DEC on platyrrhini null model"
tr = read.tree(phy)
results_object = resDECnull
scriptdir ="./R/Historical_Biogeography/a_scripts/"
# States
res2 = plot_BioGeoBEARS_results(results_object, addl_params=list("j"), plotwhat="text", label.offset=0.2, tipcex=0.35, statecex=0.25, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.2, tipcex=0.35, statecex=0.10, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#-------------------------------------------------------------------------------

#### DEC+J null model analysis ---------------------------------------------------
DECJ_platyrrhini_null = define_BioGeoBEARS_run()
DECJ_platyrrhini_null$trfn = phy
DECJ_platyrrhini_null$geogfn = table_geo
DECJ_platyrrhini_null$max_range_size = max_range_size
DECJ_platyrrhini_null$min_branchlength = 0.000001   
DECJ_platyrrhini_null$include_null_range = TRUE   
DECJ_platyrrhini_null$on_NaN_error = -1e50  
DECJ_platyrrhini_null$speedup = TRUE         
DECJ_platyrrhini_null$use_optimx = TRUE   
DECJ_platyrrhini_null$num_cores_to_use = 8
DECJ_platyrrhini_null$force_sparse = FALSE  
DECJ_platyrrhini_null$return_condlikes_table = TRUE
DECJ_platyrrhini_null$calc_TTL_loglike_from_condlikes_table = TRUE
DECJ_platyrrhini_null$calc_ancprobs = TRUE    

# Set up DEC+J model
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
DECJ_platyrrhini_null$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DECJ_platyrrhini_null$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DECJ_platyrrhini_null$BioGeoBEARS_model_object@params_table["e","init"] = estart
DECJ_platyrrhini_null$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
DECJ_platyrrhini_null$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DECJ_platyrrhini_null$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DECJ_platyrrhini_null$BioGeoBEARS_model_object@params_table["j","est"] = jstart

DECJ_platyrrhini_null = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=DECJ_platyrrhini_null)
check_BioGeoBEARS_run(DECJ_platyrrhini_null)

resfn = "./Results/Historical_biogeography/DECJ/Null/Platyrrhini_DECJ_null_v1.Rdata"
runslow = F
if (runslow){
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  res = bears_optim_run(DECJ_platyrrhini_null)
  res    
  save(res, file=resfn)
  resDECjNull = res
} else {
  load(resfn)
  resDECjNull = res
}

#######################################################
#### Plot ancestral states - DECJ null Model -----------------------------------
#######################################################
pdffn = "./Results/Historical_biogeography/DECJ/Null/Platyrrhini_DEC_null_model.pdf"
pdf(pdffn, height=15, width=7)
analysis_titletxt ="BioGeoBEARS DEC+J on platyrrhini null model"
tr = read.tree(phy)
results_object = resDECjNull
scriptdir ="./R/Historical_Biogeography/a_scripts/"
# States
res2 = plot_BioGeoBEARS_results(results_object, addl_params=list("j"), plotwhat="text", label.offset=0.2, tipcex=0.35, statecex=0.25, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.2, tipcex=0.35, statecex=0.10, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

## References ------------------------------------------------------------------
# Matzke, N.J. (2013) BioGeoBEARS: BioGeography with Bayesian (and likelihood) evolutionary analysis in R Scripts. R package, version 0.2, 1, 2013.
