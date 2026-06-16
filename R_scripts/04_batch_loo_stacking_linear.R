#-------------------------------------------------------------------------------
#
#   Stan Leave-one-out cross validation and stacking weights
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: Feb 2, 2026

# DESCRIPTION: This script loads in stanfit objects created in Stan_batch_run
# and performs Leave-one-out cross validation and estimates stacking weights. 
# Uses model stacking weights to create model stacked predicted length- and 
# growth-at-age curves at the population and sampling event-level, stacked 
# growth parameters, and isnt. growth and hydrology predictions. Additionally,
# length- and growth-at-age curves for each candidate model are created and R2s
# are estimated. This script completes all of this in batch for all species in
# analysis.


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load packages
library(dplyr)
library(tidyr)
library(loo)

# directories
input_dir <- "input_data"
fun_dir <-"functions"
out_dir <- "stan_outputs/model_out_linear"
plot_dir <- "stan_outputs/plotting_info"
export_dir <- "loo_outputs"

# Load in custom functions
source(file.path(fun_dir,"stan_loo_batch_functions.R"))
source(file.path(fun_dir,"growth_prediction_functions.R"))

# Load in data
fish_df <- readRDS(file.path(input_dir,"fsage_cleaned_2026-02-26.rds"))
sample_bridge <- readRDS(file.path(plot_dir,"fsgwh_sampleid_bridge_2026-04-23.rds"))
pred_lables <- readRDS(file.path(plot_dir,"fsgwh_pred_labels_2026-04-23.rds"))
mean_lengths <- readRDS(file.path(plot_dir,"fsgwh_mean_lengths_2026-04-23.rds"))
n.cores <- 6
stack.iter <- 10000

# Load in model outputs  -------------------------------------------------------

# Species-specific directories
sp <- list.files(out_dir)
sp_dir <- sapply(sp,function(x) file.path(out_dir,x))
names(sp_dir) <- sp

# Species specific model outputs
sp_out <- lapply(sp_dir,list.files)
sp_out <-lapply(sp_out,function(x) x[substr(x,4,6) == "111"])

# are there 3 models in each?
sapply(sp_out,function(x) length(unique(substr(x,1,6))))


# loo and model stacking  ------------------------------------------------------

# Estimate loo values for each species and model
sp_loo <- purrr::map2(sp_out,sp_dir,loo_batch,n.cores)

# Ensure that all models have Pareto' K > 0.7
sp_loo_diag <- lapply(sp_loo, loo_diag,"ESS")

# Compare loo values among models
sp_loo_compare <- lapply(sp_loo,loo_compare)

# Estimate model stacking weights
sp_stack_wt <- lapply(sp_loo,stack_format,cores = n.cores)


# Prepare prediction input data  -----------------------------------------------
# Creates sampling-event specific input data (i.e. ages 0-max.age) to be used
# to create model stacked and candidate model prediction curves

# How many sampling events exist per species?
group_size_list <- lapply(sp, function(sp) {
  g <- fish_df %>% 
    group_by(species) %>% 
    summarise(n_sample = n_distinct(wateryear,region,site)) %>% 
    filter(species == sp)
  c(g$n_sample)
})
names(group_size_list) <- sp

# Data frames containing input ages for sampling-event level predictions. 
input_bridge <- purrr::map2(
  group_size_list,names(group_size_list),
  grouping_predDF,
  min.pred=0,
  max.pred=360,
  group.id = c("")
  )
input_bridge <- purrr::transpose(input_bridge)

# Extract sampling event and age inputs
input_list <- input_bridge$prediction

# Create population input using one group idea
input_list_mu <-lapply(input_list, function(x) x %>% filter(group_id ==1))

# Extract mean age used for inst. growth predictions
mean_length_input <- lapply(
  mean_lengths, 
  function(x) data.frame(group_id=1, input =x)
  )
mean_length_input <- mean_length_input[order(names(mean_length_input))]

# DF to link group id back to organizational data
curve_id_bridge <- bind_rows(input_bridge$id_bridge)


# Create prediction arrays  ----------------------------------------------------
# Batch run functions designed to extract parameters for each candidate model
# to create prediction curves for growth and length at a range of input ages.
# First, curves are created for each candidate model, then stacking weights 
# are used to created model stacked prediction curves.

### Individual model predictions  ###

# Length and instantaneous growth at age predictions - population level
ind_mu_curve_list <- Map(function(x,d,z) 
  curve_predictR(
    stack.df = x,
    mod.dir = d,
    group.id="mu_",
    n.sim = stack.iter,
    type = "prediction",
    input.df = z,
    input.var = "age",
    output.var = c("length","growth"),
    sum.fun="median",
    parallel = T,
    mc.cores = n.cores),
  sp_stack_wt,
  sp_dir,
  input_list)

# Inst. growth at mean age
ind_mean_growth_list <- Map(function(x,d,z) 
  curve_predictR(
    stack.df = x,
    mod.dir = d,
    group.id="mu_",
    n.sim = stack.iter,
    type = "prediction",
    input.df = z,
    input.var = "length",
    output.var = c("growth"),
    sum.fun ="median",
    parallel = T,
    mc.cores = n.cores),
  sp_stack_wt,
  sp_dir,
  mean_length_input)

# Estimate R2 of candidate curves
r2_list <- lapply(1:length(sp_stack_wt),
  function(i) len_R2(
    stack.df = sp_stack_wt[[i]],
    mod.dir = sp_dir[[i]],
    group.id="",
    n.sim = stack.iter,
    sp = names(sp_stack_wt)[i],
    input.df = input_list[[i]],
    sum.fun="median",
    parallel = T,
    mc.cores = n.cores)
  )
names(r2_list) <- names(sp_stack_wt)


### Model stacked predictions  ###

# Length and instantaneous growth at age predictions - sampling-event level
curve_list <- Map(function(x,y,z) 
  growth_stackR(
    stack.df = x,
    mod.dir = y,
    group.id="",
    sim = stack.iter,
    type = "prediction",
    input.df = z,
    input.var = "age",
    output.var = c("length","growth"),
    parallel = T,
    mc.cores = n.cores,
    sum.fun="median"),
  sp_stack_wt,
  sp_dir,
  input_list)

# Length and instantaneous growth at age predictions - population level
mu_curve_list <- Map(
  function(x,y,z) growth_stackR(
    stack.df = x,
    mod.dir = y,
    group.id="mu_",
    sim = stack.iter,
    type = "prediction",
    input.df = z,
    input.var = "age",
    output.var = c("length","growth"),
    parallel = T,
    mc.cores = n.cores,
    sum.fun="median"),
  sp_stack_wt,
  sp_dir,
  input_list_mu)

# Inst. growth at mean age
mean_growth_list <- Map(
  function(x,y,z) growth_stackR(
    stack.df = x,
    mod.dir = y,
    group.id="mu_",
    sim = stack.iter,
    type = "prediction",
    input.df = z,
    input.var = "length",
    output.var = c("growth"),
    parallel = T,
    mc.cores = n.cores,
    sum.fun="median"),
  sp_stack_wt,
  sp_dir,
  mean_length_input)

# Population means of parameters
param_list <- Map(function(x,y) 
  growth_stackR(
    stack.df = x,
    mod.dir = y,
    group.id="mu_",
    sim = stack.iter,
    type = "parameter",
    truncate.inf = F,
    sum.fun="median"),
  sp_stack_wt,
  sp_dir)

# Parameter and inst. growth linear predictions
pred_list <- purrr::map2(
  sp_stack_wt,sp_dir,
  linear_pred_stackR, 
  sim = stack.iter,
  sum.fun = "median"
  )


# Format predictions into data frames  -----------------------------------------

# Add species names to data frames
pred_list <- purrr::map2(
  pred_list,names(pred_list),
  function(x,y) x %>% mutate(species =y)
  )
curve_list <- purrr::map2(
  curve_list,names(curve_list), 
  function(x,y) x %>% mutate(species =y)
  )
mu_curve_list <- purrr::map2(
  mu_curve_list,names(mu_curve_list), 
  function(x,y) x %>% mutate(species =y)
  )
ind_mu_curve_list <- purrr::map2(
  ind_mu_curve_list,names(ind_mu_curve_list),
  function(x,y) x %>% mutate(species =y)
  )
r2_list <- purrr::map2(
  r2_list,names(r2_list), 
  function(x,y) x %>% mutate(species =y)
  )
param_list <- purrr::map2(
  param_list,names(param_list), 
  function(x,y) x %>% mutate(species =y)
  )
mean_growth_list <- purrr::map2(
  mean_growth_list,names(mean_growth_list), 
  function(x,y) x %>% mutate(species =y)
  )
ind_mean_growth_list <- purrr::map2(
  ind_mean_growth_list,names(ind_mean_growth_list), 
  function(x,y) x %>% mutate(species =y)
  )

# Create one dataframe for all species
pred_df <- bind_rows(pred_list)
curve_df <- bind_rows(curve_list)
mu_curve_df <- bind_rows(mu_curve_list)
ind_mu_curve_df <- bind_rows(ind_mu_curve_list)
r2_df <- bind_rows(r2_list)
param_df <- bind_rows(param_list)
mean_growth_df <- bind_rows(mean_growth_list)
ind_mean_growth_df <- bind_rows(ind_mean_growth_list)


# Link predictions to labels  -------------------------------------------------- 
# Return sampling event labels (e.g., site and year) back to sampling-event
# level growth curves and PC values to the growth at hydrology predictions.

### Link PC values to linear predictions  ###

# Format prediction input labels to merge
pred_lables_id <- lapply(pred_lables,function(x) {
  df <- as.data.frame(x)
  df[,"pred_id"] <- 1:nrow(df)
  df
})

# add species names
pred_lables_id <- purrr::map2(
  pred_lables_id,names(pred_lables_id), 
  function(x,y) x %>% mutate(species =y)
  )

# Merge into one df
pred_lables_df <- bind_rows(pred_lables_id)

# Link prediction input labels to predictions
pred_bridged <- pred_df %>% 
  left_join(pred_lables_df) %>% 
  select(-pred_id)

### link sampling event info to curves  ###

# Add hydro and sampling info to bridge df
bridge_df <- curve_id_bridge %>% 
  left_join(sample_bridge, by = join_by(species,sample_id)) %>% 
  mutate(group_name = case_when(
    #!is.na(mu) ~ "population",
    !is.na(sample_id) ~ paste(region,site,wateryear)
  )) %>% 
  select(species,group_id,sample_id,group_name)

# Link to predictions
curve_bridged <-curve_df %>% 
  left_join(bridge_df,by = join_by(species,group_id))%>%
  select(-group_id)


# Export  ----------------------------------------------------------------------

# For summary stats, export loo, stacking weights, growth rates, r2, and parameters
saveRDS(sp_loo_compare,file.path(export_dir,paste0("loo_out_",Sys.Date(),".rds")))
saveRDS(sp_stack_wt,file.path(export_dir,paste0("stack_wt_out_",Sys.Date(),".rds")))
saveRDS(mean_growth_df, file.path(export_dir,paste0("stacked_meand_growth_predictions_",Sys.Date(),".rds")))
saveRDS(ind_mean_growth_df, file.path(export_dir,paste0("ind_mean_growth_predictions_",Sys.Date(),".rds")))
saveRDS(r2_df, file.path(export_dir,paste0("ind_model_r2_",Sys.Date(),".rds")))
saveRDS(param_df, file.path(export_dir,paste0("stacked_mu_parameters_",Sys.Date(),".rds")))

# For plots, export length and growth-at-age predictions and bridge tables
saveRDS(curve_bridged, file.path(export_dir,paste0("stacked_curves_",Sys.Date(),".rds")))
saveRDS(mu_curve_df, file.path(export_dir,paste0("stacked_mu_curves_",Sys.Date(),".rds")))
saveRDS(ind_mu_curve_df, file.path(export_dir,paste0("ind_mu_curves_",Sys.Date(),".rds")))
saveRDS(pred_bridged, file.path(export_dir,paste0("stacked_growth_predictions_",Sys.Date(),".rds")))

