library(ape)
library(optimx)
library(GenSA)
library(rexpokit)
library(cladoRcpp)
library(snow)
library(parallel)
library(BioGeoBEARS)
library(readxl)
library(tidyverse)

## Set up geography matrix -----------------------------------------------------
# Open Presence-absence matrix
presabs <- read_xlsx("../Tambouille_Lucas/BGB/Primates_geogr.xlsx")
presabs[,2:12] <- apply(X = presabs[,2:12], FUN = as.numeric, MARGIN = c(1,2))
# Remove rogue
RogueTaxa <- c("Acrecebus", "Arsinoea", "Branisella", "Catopithecus", "Chilecebus", "Panamacebus", "Parvimico",
               "Proteopithecus", "Siamopithecus", "Tremacebus")
presabs <- presabs %>% filter(Taxon %in% RogueTaxa == FALSE)
# Extract key features
n_taxa <- nrow(presabs)
n_areas <- ncol(presabs) - 1 # we don't want taxon names to be counted as areas
area_vect <- paste0("(", colnames(presabs)[2])
for(area in colnames(presabs)[3:ncol(presabs)]){
  area_vect <- paste(area_vect, area, sep = " ")
}
area_vect <- paste0(area_vect, ")")
# Convert matrix to phylip file
phylip <- paste0(nrow(presabs), "\t", n_areas, " ", area_vect, "\n")
for(i in 1:nrow(presabs)){
  # Coerce the range
  tax_range <- presabs[i, 2:12]
  coerced_range <- paste(tax_range, collapse = "")
  # Extend the phylip
  phylip <- paste0(phylip, presabs$Taxon[i], "\t", coerced_range, "\n")
}
# Save phylip
geofn <- "../Tambouille_Lucas/BGB/Primates_geogr.data"
write(phylip, "../Tambouille_Lucas/BGB/Primates_geogr.data")

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)



## Open tree
tree <- read.tree("../Tambouille_Lucas/BGB/Anthropoid_consensus_tree-10rogueless_NoHuman.nwck")

#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
tree_path <- "../Tambouille_Lucas/BGB/Anthropoid_consensus_tree-10rogueless_NoHuman.nwck"
BioGeoBEARS_run_object$trfn = tree_path

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = "../Tambouille_Lucas/BGB/Primates_geogr.data"

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = 3
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();

# n threads
BioGeoBEARS_run_object$num_cores_to_use = 10

BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "Platy_first_biogeo.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res

}

resDEC <- load("./Platy_first_biogeo.Rdata")

#######################################################
# PDF plots
#######################################################
pdffn = "Platy_unconstrained_M0.pdf"
pdf(pdffn, height=6, width=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="M0 Platy (unconstrained)"

# Setup
results_object = resDEC

scriptdir <- "../Tambouille_Lucas/BGB/"
# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, include_null_range=TRUE, cornercoords_loc = scriptdir, tr=tree, tipranges = tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
