################################################################################
# Name: 10-BioGeoBEARS_DEC_Null.R
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

#### DEC null model analysis ---------------------------------------------------
DEC_platyrrhini_null = define_BioGeoBEARS_run()
DEC_platyrrhini_null$trfn = phy
DEC_platyrrhini_null$geogfn = table_geo

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=table_geo) ## Get the maximum range size observed
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))
DEC_platyrrhini_null$max_range_size = max_range_size

areas = getareas_from_tipranges_object(tipranges) ## Get area (letter)
states_area_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=T)
DEC_platyrrhini_null$states_list = states_area_list 

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

if (sp == "all") {
  resfnN = paste("./Results/Historical_biogeography/",All,"DEC/Null_model/Platyrrhini_DECN_allprimate.Rdata",sep="")
} else if (sp == "Haplorrhini") {
  resfnN = paste("./Results/Historical_biogeography/",Hap,"DEC/Null_model/Platyrrhini_DECN_Haplo.Rdata",sep="")
} else if (sp == "Anthropoidea") {
  resfnN = paste("./Results/Historical_biogeography/",Anth,"DEC/Null_model/Platyrrhini_DECN_Anth.Rdata",sep="")
} 

runslow = T
if (runslow){
  resN= bears_optim_run(DEC_platyrrhini_null)
  resN    
  save(resN, file=resfnN)
  resDECnull = resN
} else {
  load(resfnN)
  resDECnull = resN
}

#######################################################
#### Plot ancestral states - DEC null Model ------------------------------------
#######################################################
if (sp == "all") {
  pdffn = paste("./Results/Historical_biogeography/",All,"DEC/Null_model/Platyrrhini_DECN_allprimate.pdf",sep="")
} else if (sp == "Haplorrhini") {
  pdffn = paste("./Results/Historical_biogeography/",Hap,"DEC/Null_model/Platyrrhini_DECN_Haplo.pdf",sep="")
} else if (sp == "Anthropoidea") {
  pdffn = paste("./Results/Historical_biogeography/",Anth,"DEC/Null_model/Platyrrhini_DECN_Anth.pdf",sep="")
} 
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

## References ------------------------------------------------------------------
# Matzke, N.J. (2013) BioGeoBEARS: BioGeography with Bayesian (and likelihood) evolutionary analysis in R Scripts. R package, version 0.2, 1, 2013.