################################################################################
# Name: 11-BioGeoBEARS_DEC+J_Null.R
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
  table_geo = paste("./Data/Historical_biogeography/constrained_model/",All,"Area_biogeography_platy.txt",sep="")
  
} else if (sp == "Haplorrhini") {
  phy <- paste("./Data/Historical_biogeography/Trees/All_extant/Haplo_cons_tree.newick")
  temp <- read.tree(phy)
  table_geo = paste("./Data/Historical_biogeography/constrained_model/",Hap,"Area_biogeography_platy.txt",sep="")
  
} else if (sp == "Anthropoidea") {
  phy <- paste("./Data/Historical_biogeography/Trees/All_extant/Anth_cons_tree.newick")
  temp <- read.tree(phy)
  table_geo = paste("./Data/Historical_biogeography/constrained_model/",Anth,"Area_biogeography_platy.txt",sep="")
} 


#### DEC+J null model analysis ---------------------------------------------------
DECJ_platyrrhini_null = define_BioGeoBEARS_run()
DECJ_platyrrhini_null$trfn = phy
DECJ_platyrrhini_null$geogfn = table_geo

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=table_geo) 
max(rowSums(dfnums_to_numeric(tipranges@df))) 
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))
DECJ_platyrrhini_null$max_range_size = max_range_size

areas = getareas_from_tipranges_object(tipranges) 
states_area_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=T)
DECJ_platyrrhini_null$states_list = states_area_list 

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
if (sp == "all") {
  load("./Results/Historical_biogeography/All_Primates/DEC/Null_model/Platyrrhini_DECN_allprimate.Rdata")
  resDECnull = resN
} else if (sp == "Haplorrhini") {
  load("./Results/Historical_biogeography/Haplorrhini/DEC/Null_model/Platyrrhini_DECN_Haplo.Rdata")  
  resDECnull = resN
} else if (sp == "Anthropoidea") {
  load("./Results/Historical_biogeography/Anthropoidea/DEC/Null_model/Platyrrhini_DECN_Anth.Rdata")  
  resDECnull = resN
} 
dstart = resDECnull$outputs@params_table["d","est"]
estart = resDECnull$outputs@params_table["e","est"]
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

if (sp == "all") {
  resfn_Nj = paste("./Results/Historical_biogeography/",All,"DECJ/Null_model/Platyrrhini_DECJN_allprimate.Rdata",sep="")
} else if (sp == "Haplorrhini") {
  resfn_Nj = paste("./Results/Historical_biogeography/",Hap,"DECJ/Null_model/Platyrrhini_DECJN_Haplo.Rdata",sep="")
} else if (sp == "Anthropoidea") {
  resfn_Nj = paste("./Results/Historical_biogeography/",Anth,"DECJ/Null_model/Platyrrhini_DECJN_Anth.Rdata",sep="")
} 

runslow = T
if (runslow){
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  resNj = bears_optim_run(DECJ_platyrrhini_null)
  resNj    
  save(resNj, file=resfn_Nj)
  resDECjNull = resNj
} else {
  load(resfn_Nj)
  resDECjNull = resNj
}

#######################################################
#### Plot ancestral states - DECJ null Model -----------------------------------
#######################################################
if (sp == "all") {
  pdffn = paste("./Results/Historical_biogeography/",All,"DECJ/Null_model/Platyrrhini_DECJN_allprimate.pdf",sep="")
} else if (sp == "Haplorrhini") {
  pdffn = paste("./Results/Historical_biogeography/",Hap,"DECJ/Null_model/Platyrrhini_DECJN_Haplo.pdf",sep="")
} else if (sp == "Anthropoidea") {
  pdffn = paste("./Results/Historical_biogeography/",Anth,"DECJ/Null_model/Platyrrhini_DECJN_Anth.pdf",sep="")
} 
pdf(pdffn, height=15, width=7)
analysis_titletxt ="BioGeoBEARS DEC+J null model on platyrrhini"
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