#packages#####
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
library(ltstR)
library(tidyverse)
library(dplyr)
library(ggtree)
library(igraph)

#BIOGEOBEARS####
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

# trfn <- "./iaapass.NEWICK"
trfn <- "./newtrfn.NEWICK"
geogfn <- "new.geog.DATA"

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

max(rowSums(dfnums_to_numeric(tipranges@df)))

max_range_size = 10

tr = read.tree(trfn)

# tr = read.tree(trfn)
# tr <- impose_min_brlen(tr, min_brlen = 0.000001)
# 
# write.tree(tr, file = "newtrfn.NEWICK")

######DEC######
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
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_default.txt" #NG always there
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_NG30.txt" #NG since 30Ma

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

runslow = FALSE
resfn = "iaa_DEC_tAD.Rdata"
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

##Run DEC+J########################################
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
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_default.txt" #NG always there
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_NG30.txt" #NG since 30Ma

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

resfn = "iaa_DEC+J_tAD.Rdata"
runslow = FALSE
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

##calculating AICs #####
restable2 = NULL
teststable2 = NULL

load("iaa_DEC_tND.Rdata")

LnL_2 = get_LnL_from_BioGeoBEARS_results_object(res)
# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=res, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)


load("iaa_DEC+J_tND.Rdata")

LnL_1 = get_LnL_from_BioGeoBEARS_results_object(res)

numparams1 = 5
numparams2 = 4

stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=res, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable2 = rbind(restable2, res2, res1)
teststable2 = rbind(teststable2, tmp_tests)

teststable = rbind(teststable, teststable2)
restable$x <- 0
restable$w <- 0
restable <- restable[,c(1,2,3,4,6,7,5)]
restable <- rbind(restable, restable2)
restable3 = rbin

row.names(restable) = c("DEC", "DEC+J", "BAYAREALIKE", "BAYAREALIKE+J","DEC_t", "DEC+J_t", 
                        "DEC_tA", "DEC+J_tA", "DEC_tN", "DEC+J_tN",
                        "DEC_tAD", "DEC+J_tAD", "DEC_tND", "DEC+J_tND")

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = 1390
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 <- restable
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

# save(restable, file="restable_v1.Rdata")
# save(restable2, file="restable_v2.Rdata")
save(teststable, file="teststable_v1.Rdata")
load(file="restable_v2.Rdata")
load(file="restable_v1.Rdata")
load(file="teststable_v1.Rdata")

#BSM####
model_name = "DEC+J_tAD"
res = resDECj

clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

BSM_inputs_fn = "BSM_DECj.Rdata"
runInputsSlow = FALSE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = FALSE
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  load(file="RES_clado_events_tables.Rdata")
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 10

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_matrix = col2rgb(c("#d73027", "#ffa600", "#ff7c43", "#f95d6a", "#d45087", "#a05195","#4575b4","#4575b4", "#2f4b7c","#003f5c"))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based, plot_null_range=TRUE)

##Setup for painting a single stochastic map#########
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=resDECj, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

par(mar = c(4,1,2,1))
plot_BioGeoBEARS_results(results_object=resDECj, analysis_titletxt="DEC+J_tAD", 
                         addl_params=list("j"), label.offset=0.8, plotwhat="pie", 
                         tipcex = 0.5, statecex = 0.2, plotsplits = FALSE,
                         cornercoords_loc=scriptdir, root.edge=TRUE, 
                         colors_list_for_states=colors_list_for_states, 
                         skiptree=FALSE, show.tip.label=FALSE)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="DEC+J_tAD_BSM", 
                         addl_params=list("j"), label.offset=0.8, plotwhat="pie", 
                         tipcex = 0.5, statecex = 0.3, plotsplits = FALSE,
                         cornercoords_loc=scriptdir, root.edge=TRUE, 
                         colors_list_for_states=colors_list_for_states, 
                         skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, 
                              colors_list_for_states=colors_list_for_states, lwd=2, lty=par("lty"), 
                              root.edge=TRUE, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="DEC+J_tAD_BSM", 
                         addl_params=list("j"), label.offset=0.8,plotwhat="pie", 
                         cornercoords_loc=scriptdir, tipcex = 0.5, statecex = 0.3, plotsplits = FALSE,
                         root.edge=TRUE, colors_list_for_states=colors_list_for_states, 
                         skiptree=TRUE, show.tip.label=TRUE)

##Summarize stochastic map tables#################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

#Counting Number of Events total ####
counts_list = count_ana_clado_eventsM(clado_events_tables, ana_events_tables, areanames, actual_names)

BSM_sum <- as.data.frame(t(counts_list$summary_counts_BSMs))
BSM_sum_sum <- BSM_sum[c(8, 6, 5, 7),]
BSM_sum_sum[4,] <- BSM_sum_sum[3,] + BSM_sum_sum[4,]
BSM_sum_sum <- BSM_sum_sum[-3,]
BSM_sum_sum$label <- c("Dispersals", "Vicariance", "Sympatry")
BSM_sum_sum$label <- factor(BSM_sum_sum$label, levels = c("Dispersals", "Sympatry", "Vicariance"))

ggplot(data=BSM_sum_sum, aes(x=label, y=means, fill = label)) +
  geom_bar(stat = "identity") +
  xlab("Event Type") +
  ylab("Mean No of Events") +
  theme_minimal()

BSM_sum <- BSM_sum[c(1,3,5,6,7),]
BSM_sum$label <- c("founder", "expansion", "subset", "vicariance", "sympatry")
BSM_sum$label <- factor(BSM_sum$label, levels = c("expansion", "founder", "sympatry", "subset", "vicariance"))
# 
eventscountsplot <- ggplot(data=BSM_sum, aes(x=label, y=means, fill = label)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("b) Event Counts") +
  xlab("") +
  ylab("mean no of events") +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

# Identifying back-colonization ####
## How many upstream vs downstream colonisation? 
total.disp <- counts_list$all_dispersals_counts_fromto_means
island.lab <- c("E", "M", "U", "L")
cont.lab <- c("A", "N", "B", "J", "S", "O")
sahul.lab <- c("A", "N")
sunda.lab <- c("B", "J", "S", "O")

total.disp.wself <- sum(total.disp) - sum(total.disp[island.lab, island.lab]) - sum(total.disp[cont.lab, cont.lab])
sum(total.disp[island.lab, cont.lab])/sum(total.disp.wself) #upstream
sum(total.disp[cont.lab, island.lab])/sum(total.disp.wself) #downstream
sum(total.disp[island.lab, island.lab]) #within island
sum(total.disp[cont.lab, cont.lab])

###aggregating... ####
total.disp.agg <- 
  cbind(rowSums(total.disp[, 1:2]),
      total.disp[,3],
      rowSums(total.disp[, 4:6]),
      rowSums(total.disp[, 7:10]))

total.disp.agg <- 
  rbind(colSums(total.disp.agg[1:2, ]), 
        total.disp.agg[3, ],
        colSums(total.disp.agg[4:6, ]),
        colSums(total.disp.agg[7:10, ]))

colnames(total.disp.agg) <- c("Sahul", "EMN", "Wallacea", "Sunda")
rownames(total.disp.agg) <- c("Sahul", "EMN", "Wallacea", "Sunda")

bar.disp <- melt(total.disp.agg)
colnames(bar.disp) <- c("from","to","disp")
bar.disp$from <- factor(bar.disp$from, levels = c("Sunda", "Wallacea", "EMN", "Sahul"))
bar.disp$to <- factor(bar.disp$to, levels = c("Sunda", "Wallacea", "EMN", "Sahul"))

bar.disp <- bar.disp[bar.disp$from != bar.disp$to, ]

ggplot(data = bar.disp, aes(fill= to, y = disp, x = from)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c( "#00b050","#48acc8","#43e9d8", "#fd696c")) +
  scale_y_continuous(limits = c(0, 200)) +
  theme_minimal()+
  xlab("Sources")+
  ylab("No. of outward dispersals") +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

## Founder Events ####
### which nodes are founder events upstream####
tr <- read.tree(trfn)

f.nodes.summ <- NULL
i=0
id.founder <- 
  function(clado_events_tables){
    for(i in 1:length(clado_events_tables)){
      f.nodes <- subset(clado_events_tables[[i]],clado_events_tables[[i]]$clado_event_type=="founder (j)", 
                        select=c(node, time_bp, clado_dispersal_from, clado_dispersal_to))
      f.nodes <- cbind(f.nodes,rep(i, nrow(f.nodes)))
      
      f.nodes <- subset( f.nodes, f.nodes$clado_dispersal_from=="M" | f.nodes$clado_dispersal_from=="E" |
                           f.nodes$clado_dispersal_from=="L" | f.nodes$clado_dispersal_from=="U") #make this only upstream founder
      
      f.nodes.summ <- rbind(f.nodes.summ, f.nodes)
      i=i+1
    }
    
    id.f.nodes <- aggregate(f.nodes.summ$`rep(i, nrow(f.nodes))`, by=list(f.nodes.summ$node), FUN="length")
    colnames(id.f.nodes) <- c("node", "freq")
    
    nodelabellist <- prt(tr, printflag = TRUE, get_tipnames=TRUE)
    
    id.f.nodes <- merge(id.f.nodes, nodelabellist, by = "node", all.x = TRUE, all.y = FALSE)
    id.f.nodes <- id.f.nodes[order(id.f.nodes$freq),]
    return(id.f.nodes)
  }

wherefounder <- id.founder(clado_events_tables)

## ID d anagenetic dispersal
d.nodes.summ <- NULL
i=0
id.ddisp <- 
  function(ana_events_tables){
    for(i in 1:length(ana_events_tables)){
      d.nodes <- subset(ana_events_tables[[i]],ana_events_tables[[i]]$event_type=="d", 
                        select=c(node, time_bp, event_txt))
      d.nodes <- separate(d.nodes, event_txt, into = c("from", "to"), sep = "->")
      d.nodes$island <- ifelse(grepl("^[MELU]+$", d.nodes$from), "yes", "no") #id island disp
      
      d.nodes <- cbind(d.nodes,rep(i, nrow(d.nodes)))
      
      d.nodes <- subset( d.nodes, d.nodes$island == "yes") #make this only upstream founder
      d.nodes <- d.nodes[!duplicated(d.nodes$node), ] #ensure we only want no of maps with the events, not adding up events
      
      d.nodes.summ <- rbind(d.nodes.summ, d.nodes)
      i=i+1
    }
    
    id.d.nodes <- aggregate(d.nodes.summ$`rep(i, nrow(d.nodes))`, by=list(d.nodes.summ$node), FUN="length")
    colnames(id.d.nodes) <- c("node", "freq")
    
    nodelabellist <- prt(tr, printflag = TRUE, get_tipnames=TRUE)
    
    id.d.nodes <- merge(id.d.nodes, nodelabellist, by = "node", all.x = TRUE, all.y = FALSE)
    id.d.nodes <- id.d.nodes[order(id.d.nodes$freq),]
    return(id.d.nodes)
  }

whereddisp <- id.ddisp(ana_events_tables)


timebp <- data.frame(node = clado_events_tables[[1]]$node, timebp = clado_events_tables[[1]]$time_bp)
### Plotting ####
tr.dt <- fortify(tr)
tr.dt <- merge(tr.dt, timebp, by = "node")
tr.dt <- merge(tr.dt, wherefounder[, c("node", "freq")], by = "node", all.x = TRUE)
names(tr.dt)[names(tr.dt) == "freq"] <- "founder"

# intervals <- c(5, 16, 25, 50)
# labels <- c("low", "medium", "high")
# tr.dt$founder <- cut(tr.dt$founder, breaks = intervals, labels = labels, include.lowest = TRUE)

tr.dt <- merge(tr.dt, whereddisp[, c("node", "freq")], by = "node", all.x = TRUE)
names(tr.dt)[names(tr.dt) == "freq"] <- "d"
# intervals <- c(4, 12, 25, 50)
# labels <- c("low", "medium", "high")
# tr.dt$d <- cut(tr.dt$d, breaks = intervals, labels = labels, include.lowest = TRUE)

# tr.dt <- tr.dt %>%
#   mutate(disp = case_when(
#     !is.na(d) & is.na(founder) ~ d,  # if only 'a' has a value
#     is.na(d) & !is.na(founder) ~ founder,  # if only 'b' has a value
#     !is.na(d) & !is.na(founder) ~ ifelse(d < founder, founder, d)  # if both 'a' and 'b' have values
#   ))

tr.dt$disp <- tr.dt$founder + tr.dt$d

intervals <- c(quantile(tr.dt$disp, probs = c(0.25, 0.5, 0.75, 1), na.rm = TRUE))
labels <- c("low", "medium", "high")
tr.dt$disp <- cut(tr.dt$disp, breaks = intervals, labels = labels, include.lowest = TRUE)
testing <- c("high", "medium", "low")

for (label in testing) {
  # Get nodes corresponding to the current category
  nodes <- tr.dt$node[tr.dt$disp == label]
  
  # Iterate through nodes of the current category
  for (node in nodes[!is.na(nodes)]) {
    # Get descendants of the current node
    descendants <- unlist(getDescendants(tr, node))
    
    # Update category for descendants
    tr.dt$disp[tr.dt$node %in% descendants] <- label
  }
}

tr.dt %>% 
  filter(timebp >= 0 & timebp <= 30) %>% 
  group_by(disp) %>% 
  summarise(count = n())

ggplot(tr.dt[!is.na(tr.dt$disp), ], aes(x = cut(timebp, breaks = seq(0,60, by =5)), fill = disp, color = disp)) +
  geom_bar(position = "stack") +
  theme_classic()+
  xlab("Time (myr)")+
  ylab("island-to-continent dispersals") +
  scale_y_continuous(expand = c(0,0), limits = c(0,700), breaks = seq(0,700, by = 100)) +
  theme( panel.grid.minor = element_blank(), legend.position = "none")

#### TREE ####
tip_source <- list(
  ACANTHIZIDAE=c('Acanthiza_cinerea', 'Pachycare_flavogriseum'),
  PARDALOTIDAE=c('Pardalotus_striatus', 'Pardalotus_rubricatus'),
  DASYORNITHIDAE=c('Dasyornis_broadbenti', 'Dasyornis_brachypterus'),
  MALURIDAE=c('Malurus_alboscapulatus', 'Amytornis_barbatus'),
  POMATOSTOMIDAE=c('Pomatostomus_superciliosus', 'Pomatostomus_ruficeps'),
  PTILONORHYNCHIDAE=c('Ailuroedus_buccoides', 'Ptilonorhynchus_violaceus'),
  CLIMACTERIDAE=c('Climacteris_picumnus','Cormobates_leucophaea'),
  ATRICHORNITHIDAE=c('Atrichornis_rufescens', 'Atrichornis_clamosus'),
  MENURIDAE=c('Menura_novaehollandiae', 'Menura_alberti'),
  PITTIDAE=c('Erythropitta_granatina','Pitta_superba'),#'Pitta_erythrogaster',
  CALYPTOMENIDAE=c('Calyptomena_viridis', 'Calyptomena_hosii'),
  EURYLAIMIDAE=c('Eurylaimus_ochromalus', 'Corydon_sumatranus'),
  CINCLOSOMATIDAE=c('Cinclosoma_punctatum','Ptilorrhoa_leucosticta','Cinclosoma_marginatum'),
  CAMPEPHAGIDAE=c('Edolisoma_tenuirostre','Pericrocotus_divaricatus', 'Pericrocotus_montanus'),
  NEOSITTIDAE=c('Daphoenositta_chrysoptera', 'Daphoenositta_miranda'),
  PSOPHODIDAE=c('Psophodes_cristatus', 'Androphobus_viridis'),
  OREOICIDAE=c('Oreoica_gutturalis', 'Ornorectes_cristatus'),
  ORIOLIDAE=c('Oriolus_chinensis','Pitohui_dichrous','Sphecotheres_hypoleucus'),
  PARAMYTHIIDAE=c('Oreocharis_arfaki', 'Paramythia_olivacea'),
  VIREONIDAE=c('Erpornis_zantholeuca', 'Pteruthius_aeralatus'),
  MACHAERIRHYNCHIDAE=c('Machaerirhynchus_nigripectus','Machaerirhynchus_flaviventer'),
  ARTAMIDAE=c('Artamus_cinereus', 'Melloria_quoyi', 'Peltops_blainvillii', 'Strepera_graculina'),
  AEGITHINIDAE=c('Aegithina_lafresnayei', 'Aegithina_viridissima'),
  VANGIDAE=c('Tephrodornis_virgatus', 'Philentoma_velata'),
  RHIPIDURIDAE=c('Rhipidura_javanica','Chaetorhynchus_papuensis','Rhipidura_tenebrosa'),
  DICRURIDAE=c('Dicrurus_aeneus', 'Dicrurus_sumatranus'),
  MONARCHIDAE=c('Hypothymis_azurea', 'Monarcha_megarhynchus', 'Arses_insularis', 'Myiagra_eichhorni'),
  PARADISAEIDAE=c('Paradisaea_minor', 'Phonygammus_keraudrenii'),
  CORCORACIDAE=c('Struthidea_cinerea', 'Corcorax_melanorhamphos'),
  MELAMPITTIDAE=c('Melampitta_lugubris', 'Megalampitta_gigantea'),
  CNEMOPHILIDAE=c('Cnemophilus_loriae', 'Loboparadisea_sericea'),
  MELANOCHARITIDAE=c('Melanocharis_versteri', 'Oedistoma_iliolophus', 'Toxorhampus_novaeguineae'),
  PETROICIDAE=c('Devioeca_papuana', 'Petroica_multicolor', 'Amalocichla_incerta'),
  STENOSTIRIDAE=c('Chelidorhynx_hypoxanthus', 'Culicicapa_ceylonensis'),
  ORTHONYCHIDAE=c('Orthonyx_temminckii', 'Orthonyx_novaeguineae'),
  ACROCEPHALIDAE=c('Acrocephalus_orientalis', 'Arundinax_aedon'),
  LOCUSTELLIDAE=c('Locustella_lanceolata', 'Megalurulus_grosvenori', 'Poodytes_gramineus'),#'Robsonius_rabori', 
  HIRUNDINIDAE=c('Hirundo_rustica', 'Delichon_lagopodum', 'Cecropis_daurica'),
  PYCNONOTIDAE=c('Pycnonotus_jocosus', 'Alophoixus_pallidus'),
  ZOSTEROPIDAE=c('Zosterops_everetti', 'Yuhina_castaniceps'),
  TIMALIIDAE=c('Timalia_pileata', 'Stachyris_grammiceps'),
  PELLORNEIDAE=c('Trichastoma_pyrrogenys', 'Turdinus_atrigularis'), #'Leonardina_woodi'
  SCOTOCERCIDAE=c('Abroscopus_albogularis','Tesia_cyaniventer', 'Tesia_olivea'),
  DICAEIDAE=c('Dicaeum_aureolimbatum', 'Prionochilus_maculatus'), #'Dicaeum_hypoleucum', 'Dicaeum_bicolor'
  NECTARINIIDAE=c('Leptocoma_sperata', 'Cinnyris_idenburgi'),
  #IRENIDAE=c('Irena_puella'),#'Irena_cyanogastra', 
  CHLOROPSEIDAE=c('Chloropsis_cochinchinensis', 'Chloropsis_sonnerati'),
  PASSERIDAE=c('Passer_domesticus', 'Passer_montanus' ),#'Hypocryptadius_cinnamomeus'
  MOTACILLIDAE=c('Motacilla_alba', 'Anthus_godlewskii'),
  FRINGILLIDAE=c('Fringilla_montifringilla', 'Pyrrhula_nipalensis'),#'Pyrrhula_leucogenis'
  MELIPHAGIDAE=c('Meliphaga_lewinii','Acanthorhynchus_tenuirostris', 'Myzomela_cruentata'),
  PACHYCEPHALIDAE=c('Colluricincla_harmonica','Pachycephala_simplex'),
  CORVIDAE=c('Corvus_coronoides', 'Cissa_chinensis', 'Platysmurus_leucopterus'),
  LANIIDAE=c('Lanius_schach', 'Lanius_cristatus'),
  PARIDAE=c('Parus_cinereus', 'Melanochlora_sultanea'), #'Pardaliparus_elegans'
  ALAUDIDAE=c('Alauda_gulgula', 'Mirafra_javanica'),
  LEIOTHRICHIDAE=c('Leiothrix_laurinae', 'Garrulax_castanotis'),
  PHYLLOSCOPIDAE=c('Phylloscopus_trivirgatus', 'Phylloscopus_trochiloides', 'Phylloscopus_fuscatus'),
  TURDIDAE=c('Cochoa_viridis', 'Zoothera_aurea'), #'Turdus_chrysolaus', 
  MUSCICAPIDAE=c('Muscicapa_sibirica', 'Larvivora_ruficeps', 'Saxicola_maurus'),
  STURNIDAE=c('Sturnia_malabarica', 'Aplonis_minor'), #'Rhabdornis_inornatus'
  SITTIDAE=c('Sitta_frontalis', 'Sitta_neglecta'),
  PLOCEIDAE=c('Ploceus_philippinus', 'Ploceus_manyar'),
  ESTRILDIDAE = c('Poephila_personata', 'Lonchura_leucogastra', 'Erythrura_trichroa'),
  EMBERIZIDAE = c('Emberiza_fucata', 'Emberiza_rutila'),#'Emberiza_sulphurata'
  CISTICOLIDAE = c('Orthotomus_sepium','Cisticola_exilis'),#'Orthotomus_castaneiceps', 
  SYLVIIDAE = c('Cholornis_unicolor', 'Paradoxornis_guttaticollis')
)

mrca_nd <- data.frame(Family = character(), MRCA_Node = integer(), stringsAsFactors = FALSE)

# Iterate through each family
for (family in names(tip_source)) {
  tip_names <- tip_source[[family]]
  
  # Find the MRCA for the current family
  mrca <- getMRCA(tr, tip_names)
  
  # Add the results to the dataframe
  mrca_nd <- rbind(mrca_nd, data.frame(Family = family, MRCA_Node = mrca, stringsAsFactors = FALSE))
}

mrca_nd <- mrca_nd[order(mrca_nd$Family), ]

mrca.nodes <- mrca_nd$MRCA_Node

ggtree(tr, layout ="fan", open.angle = 21.9) +
  geom_tiplab(geom = "circle", node = mrca.nodes, size = 2)

ggtree(tr, layout ="fan", open.angle = 21.9) +
  geom_tree(aes(color = disp), data = tr.dt) +
  theme(legend.position = "none") +
  geom_cladelab(node = MRCA(tr, "Acanthiza_cinerea", "Pachycare_flavogriseum"),
                label = "ACANTHIZIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Pardalotus_striatus", "Pardalotus_rubricatus"),
                label = "PARDALOTIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Dasyornis_broadbenti", "Dasyornis_brachypterus"),
                label = "DASYORNITHIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Malurus_alboscapulatus", "Amytornis_barbatus"),
                label = "MALURIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Pomatostomus_superciliosus", "Pomatostomus_ruficeps"),
                label = "POMATOSTOMIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Ailuroedus_buccoides", "Ptilonorhynchus_violaceus"),
                label = "PTILONORHYNCHIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Climacteris_picumnus", "Cormobates_leucophaea"),
                label = "CLIMACTERIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Atrichornis_rufescens", "Atrichornis_clamosus"),
                label = "ATRICHORNITHIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Menura_novaehollandiae", "Menura_alberti"),
                label = "MENURIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Erythropitta_granatina", "Pitta_superba"),
                label = "PITTIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Calyptomena_viridis", "Calyptomena_hosii"),
                label = "CALYPTOMENIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Eurylaimus_ochromalus", "Corydon_sumatranus"),
                label = "EURYLAIMIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Cinclosoma_punctatum", "Cinclosoma_marginatum"),
                label = "CINCLOSOMATIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Edolisoma_tenuirostre', 'Pericrocotus_divaricatus'),
                label = "CAMPEPHAGIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Daphoenositta_chrysoptera", "Daphoenositta_miranda"),
                label = "NEOSITTIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Psophodes_cristatus", "Androphobus_viridis"),
                label = "PSOPHODIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Oreoica_gutturalis", "Ornorectes_cristatus"),
                label = "OREOICIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Oriolus_chinensis', 'Pitohui_dichrous'),
                label = "ORIOLIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Oreocharis_arfaki", "Paramythia_olivacea"),
                label = "PARAMYTHIIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Erpornis_zantholeuca", "Pteruthius_aeralatus"),
                label = "VIREONIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Machaerirhynchus_nigripectus", "Machaerirhynchus_flaviventer"),
                label = "MACHAERIRHYNCHIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Artamus_cinereus', 'Strepera_graculina'),
                label = "ARTAMIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Aegithina_lafresnayei", "Aegithina_viridissima"),
                label = "AEGITHINIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Tephrodornis_virgatus", "Philentoma_velata"),
                label = "VANGIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Rhipidura_javanica", "Rhipidura_tenebrosa"),
                label = "RHIPIDURIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Dicrurus_aeneus", "Dicrurus_sumatranus"),
                label = "DICRURIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Hypothymis_azurea', 'Monarcha_megarhynchus'),
                label = "MONARCHIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Paradisaea_minor", "Phonygammus_keraudrenii"),
                label = "PARADISAEIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Struthidea_cinerea", "Corcorax_melanorhamphos"),
                label = "CORCORACIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Melampitta_lugubris", "Megalampitta_gigantea"),
                label = "MELAMPITTIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Cnemophilus_loriae", "Loboparadisea_sericea"),
                label = "CNEMOPHILIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Melanocharis_versteri", "Toxorhampus_novaeguineae"),
                label = "MELANOCHARITIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Devioeca_papuana", "Petroica_multicolor"),
                label = "PETROICIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Chelidorhynx_hypoxanthus", "Culicicapa_ceylonensis"),
                label = "STENOSTIRIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Orthonyx_temminckii", "Orthonyx_novaeguineae"),
                label = "ORTHONYCHIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Acrocephalus_orientalis", "Arundinax_aedon"),
                label = "ACROCEPHALIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Locustella_lanceolata", 'Megalurulus_grosvenori'),
                label = "LOCUSTELLIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Hirundo_rustica", "Cecropis_daurica"),
                label = "HIRUNDINIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Pycnonotus_jocosus", "Alophoixus_pallidus"),
                label = "PYCNONOTIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Zosterops_everetti", "Yuhina_castaniceps"),
                label = "ZOSTEROPIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Timalia_pileata", "Stachyris_grammiceps"),
                label = "TIMALIIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Trichastoma_pyrrogenys", "Turdinus_atrigularis"),
                label = "PELLORNEIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Abroscopus_albogularis", "Tesia_cyaniventer"),
                label = "SCOTOCERCIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Dicaeum_aureolimbatum", "Prionochilus_maculatus"),
                label = "DICAEIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Leptocoma_sperata", "Cinnyris_idenburgi"),
                label = "NECTARINIIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Chloropsis_cochinchinensis", "Chloropsis_sonnerati"),
                label = "CHLOROPSEIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Passer_domesticus", "Passer_montanus"),
                label = "PASSERIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Motacilla_alba", "Anthus_godlewskii"),
                label = "MOTACILLIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Fringilla_montifringilla", "Pyrrhula_nipalensis"),
                label = "FRINGILLIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Meliphaga_lewinii', 'Acanthorhynchus_tenuirostris'),
                label = "MELIPHAGIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Colluricincla_harmonica", "Pachycephala_simplex"),
                label = "PACHYCEPHALIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Corvus_coronoides', 'Cissa_chinensis', 'Platysmurus_leucopterus'),
                label = "CORVIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Lanius_schach", "Lanius_cristatus"),
                label = "LANIIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Parus_cinereus", "Melanochlora_sultanea"),
                label = "PARIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Alauda_gulgula", "Mirafra_javanica"),
                label = "ALAUDIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Leiothrix_laurinae", "Garrulax_castanotis"),
                label = "LEIOTHRICHIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Phylloscopus_trivirgatus', 'Phylloscopus_trochiloides', 'Phylloscopus_fuscatus'),
                label = "PHYLLOSCOPIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Cochoa_viridis", "Zoothera_aurea"),
                label = "TURDIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Muscicapa_sibirica', 'Larvivora_ruficeps', 'Saxicola_maurus'),
                label = "MUSCICAPIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Sturnia_malabarica", "Aplonis_minor"),
                label = "STURNIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Sitta_frontalis", "Sitta_neglecta"),
                label = "SITTIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Ploceus_philippinus", "Ploceus_manyar"),
                label = "PLOCEIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, 'Poephila_personata', 'Lonchura_leucogastra', 'Erythrura_trichroa'),
                label = "ESTRILDIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Emberiza_fucata", "Emberiza_rutila"),
                label = "EMBERIZIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Orthotomus_sepium", "Cisticola_exilis"),
                label = "CISTICOLIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0) +
  geom_cladelab(node = MRCA(tr, "Cholornis_unicolor", "Paradoxornis_guttaticollis"),
                label = "SYLVIIDAE\n", barsize = 3, offset = 0.5, align = T, angle = 0)
  

#Lineages through time plots ####
# create a list of bLTT curve data as we did above for each single stochastic map
CET <- clado_events_tables
AET <- ana_events_tables

CAET <- correctNodeHeights(CET,AET)
CET <- CAET[[1]]
AET <- CAET[[2]]

LTST_dt_list <- getLTSTFromStochasticMaps(tr, CET, AET, tipranges)

window_dt_list <- lapply(LTST_dt_list, slidingWindowForStochasticMaps, window= 5)

# finally to plot the range of vlaues from across stochastic maps we want ot first summarise them
# here we use the upper (97.5%) and lower (2.5%) quantiles
# melt = T presnet the output in a way that is useful for plotting with ggplot2

LTST_dt_quantiles <- timeRangeAcrossStochasticMaps(tr, window_dt_list, melt=T)

#window_dt_list5 <- lapply(LTST_dt_list, slidingWindowForStochasticMaps, window= 5)
#LTST_dt_quantiles5 <- timeRangeAcrossStochasticMaps(tr, window_dt_list5, melt=T)
#write.csv(LTST_dt_quantiles5, "lineages5myaUWMA.csv")



# and again on the log scale
# first however have to get rid of values which are -inf which results from taking log(0)
LTST_dt_quantiles$logU <- log(LTST_dt_quantiles$upperQ)
LTST_dt_quantiles$logL <- log(LTST_dt_quantiles$lowerQ)
LTST_dt_quantiles$logU[which(LTST_dt_quantiles$logU == log(0))] <- 0
LTST_dt_quantiles$logL[which(LTST_dt_quantiles$logL == log(0))] <- 0
#

LTST_dt_quantiles$region <- factor(LTST_dt_quantiles$region,
                                   levels = c("O", "S", "B", "J", "L", "U", "M", "E", "N", "A"),
                                   labels = c("SEA","SUM", "BOR", "JAV","LSU", "SUL","MAU", "EMN","NWG", "AUS"))

LTST_dt_quantiles2 <- LTST_dt_quantiles %>% 
  mutate(group = case_when(
    region %in% c("SEA", "SUM", "BOR", "JAV") ~ "Sunda",
    region %in% c("LSU", "SUL", "MAU","EMN") ~ "Islands", 
    region %in% c("NWG", "AUS") ~ "Sahul", 
    TRUE ~ NA_character_))

# plot range of data from stochastic mapping
ggplot(aes(x=time, ymin=lowerQ, ymax=upperQ, colour=region, fill=region), data=LTST_dt_quantiles2)+ 
  geom_ribbon(alpha = 5/10)+
  theme_classic()+
  xlab("Time (myr)")+
  ylab("Species Diversity")+
  scale_y_continuous(expand = c(0,0)) +
  scale_x_reverse() +
  scale_fill_manual(values=rev(hex)) +
  scale_color_manual(values=rev(hex)) +
  facet_wrap(~group)

lttlog <- ggplot(aes(x=time, ymin=logL, ymax=logU, colour=region, fill=region), data=LTST_dt_quantiles2)+ 
  geom_ribbon(alpha = 5/10)+
  ggtitle("c) Species Diversity")+
  facet_wrap(~group, ncol = 3) +
  scale_fill_manual(values=rev(hex)) +
  scale_color_manual(values=rev(hex)) +
  theme_classic()+
  xlab("Time (myr)")+
  ylab("mean and quantile (log)")+
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
  scale_x_reverse() +
  theme(legend.position = "none", axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
