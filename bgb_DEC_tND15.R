library(optimx)
library(FD)
library(snow)
library(parallel)
library(GenSA)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(phytools)
#library(ltstR)
library(tidyverse)

#BIOGEOBEARS
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

# trfn <- "./iaapass.NEWICK"
trfn <- "newtrfn.NEWICK"
geogfn <- "new.geog.DATA"

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

max(rowSums(dfnums_to_numeric(tipranges@df)))

max_range_size = 10

# tr = read.tree(trfn)
# tr <- impose_min_brlen(tr, min_brlen = 0.000001)
# 
# write.tree(tr, file = "newtrfn.NEWICK")

#####
#DEC#
#####
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

BioGeoBEARS_run_object$trfn = trfn

BioGeoBEARS_run_object$geogfn = geogfn

BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
BioGeoBEARS_run_object$num_cores_to_use = 20

BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

#time stratification
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"

#areas allowed, change timing of NG
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_default.txt" #NG always there
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_NG15.txt" #NG since 15Ma

#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"

#paleodistance included
BioGeoBEARS_run_object$distsfn = "distances_calc.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multiplier.txt" #include this to consider 'area touching'

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

#Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

BioGeoBEARS_run_object

BioGeoBEARS_run_object$BioGeoBEARS_model_object

# Add distance-dependent dispersal (+x)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 300

# Add manual dispersal multiplier (+w)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","min"] = -1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","max"] = 10

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-12
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 10

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-12
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 10

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "iaa_DEC_tND15.Rdata"
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

#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();

BioGeoBEARS_run_object$num_cores_to_use = 20
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

#time stratification
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"

#areas allowed, change timing of NG
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_default.txt" #NG always there
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_NG15.txt" #NG since 15Ma

#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"

#paleodistance included
BioGeoBEARS_run_object$distsfn = "distances_calc.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multiplier.txt" #include this to consider 'area touching'

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

#Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 300

# Add manual dispersal multiplier (+w)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","min"] = -1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","max"] = 10

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-12
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 10

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-12
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 10

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "iaa_DEC+J_tND15.Rdata"
runslow = TRUE
if (runslow)
{
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}