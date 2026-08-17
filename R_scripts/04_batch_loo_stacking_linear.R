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
out_dir <- "stan_outputs/model_out"
plot_dir <- "stan_outputs/plotting_info"
export_dir <- "loo_outputs"

# Load in custom functions
# source(file.path(fun_dir,"stan_loo_batch_functions.R"))
# source(file.path(fun_dir,"growth_prediction_functions.R"))
devtools::load_all("~/Documents/work/R packages/growthstack")

# Load in data
fish_df <- 
  readRDS(file.path(input_dir,"fsage_cleaned_2026-06-18.rds"))
sample_bridge <- 
  readRDS(file.path(plot_dir,"fsgwh_sampleid_bridge_2026-06-16.rds"))
pred_lables <- 
  readRDS(file.path(plot_dir,"fsgwh_pred_labels_2026-06-16.rds"))
mean_lengths <- 
  readRDS(file.path(plot_dir,"fsgwh_mean_lengths_2026-06-16.rds"))
n.cores <- 6
stack.iter <- 10000

# Load in model outputs  -------------------------------------------------------

# Species-specific directories
sp <- list.files(out_dir)
sp_dir <- sapply(sp,function(x) file.path(out_dir,x))
names(sp_dir) <- sp

# Species specific model outputs. Select only linear or random models
sp_out <- lapply(sp_dir,list.files,pattern = "linear|random") 

# are there 3 models in each?
sapply(sp_out,n_distinct)

# Format mean lengths
mean_lengths <- mean_lengths[order(names(mean_lengths))]


# loo and model stacking  ------------------------------------------------------

# Estimate loo values for each species and model
sp_loo <- purrr::map2(sp_out,sp_dir,loo_batch,n.cores)

# Ensure that all models have Pareto' K > 0.7
sp_loo_diag <- lapply(sp_loo, loo_diag,"ESS")

# Compare loo values among models
sp_loo_compare <- lapply(sp_loo,loo_compare)

# Estimate model stacking weights
sp_stack_wt <- lapply(sp_loo,stack_format,cores = n.cores)


# Create prediction arrays  ----------------------------------------------------
# Batch run functions designed to extract parameters for each candidate model
# to create prediction curves for growth and length at a range of input ages.
# First, curves are created for each candidate model, then stacking weights 
# are used to created model stacked prediction curves.

pred_input <- 0:360

### Individual model predictions  ###

# Length and instantaneous growth at age predictions - population level
ind_mu_curve_list <- purrr::map2(
  sp_stack_wt,
  sp_dir,
  stack_predict,
    type = "prediction",
    group.id="mu",
    sim = stack.iter,
    stack=F,
    pred.input = pred_input,
    input.var = "age",
    output.var = c("length","growth"),
    summarize=T,
    sum.fun="median",
    parallel = T,
    mc.cores = n.cores
)

# Inst. growth at mean age
ind_mean_growth_list <- Map(function(x,y,z) 
  stack_predict(
    stack.df = x,
    mod.dir = y,
    type = "prediction",
    group.id="mu",
    sim = stack.iter,
    stack=F,
    pred.input = z,
    input.var = "length",
    output.var = "growth",
    sum.fun="median",
    parallel = T,
    mc.cores = n.cores
  ),
  sp_stack_wt,
  sp_dir,
  mean_lengths
)


### Model stacked predictions  ###

# Length and instantaneous growth at age predictions - sampling-event level
site_curve_list <- purrr::map2(
  sp_stack_wt,
  sp_dir,
  stack_predict,
    group.id="site",
    sim = stack.iter,
    stack=T,
    type = "prediction",
    pred.input = pred_input,
    input.var = "age",
    output.var = c("length","growth"),
    summarize=T,
    sum.fun="median",
    parallel = T,
    mc.cores = n.cores
)

# Length and instantaneous growth at age predictions - population level
mu_curve_list <- purrr::map2(
  sp_stack_wt,
  sp_dir,
  stack_predict,
    group.id = "mu",
    sim = stack.iter,
    stack = T,
    type = "prediction",
    pred.input = pred_input,
    input.var = "age",
    output.var = c("length","growth"),
    summarize=T,
    sum.fun = "median",
    parallel = T,
    mc.cores = n.cores
)

# Inst. growth at mean age
mean_growth_list <- Map(function(x,y,z) 
  stack_predict(
    stack.df = x,
    mod.dir = y,
    type = "prediction",
    group.id="mu",
    sim = stack.iter,
    stack=T,
    pred.input = z,
    input.var = "length",
    output.var = "growth",
    summarize=T,
    sum.fun="median",
    parallel = T,
    mc.cores = n.cores
  ),
  sp_stack_wt,
  sp_dir,
  mean_lengths
)

# Population means of parameters
param_list <- purrr::map2(
  sp_stack_wt,
  sp_dir,
  stack_predict,
    group.id = "mu",
    sim = stack.iter,
    stack = T,
    type = "parameter",
    truncate.inf = F,
    summarize=T,
    sum.fun = "median"
)

# Parameter and inst. growth linear predictions
pred_list <- purrr::map2(
  sp_stack_wt[names(sp_stack_wt) !="JORFLO"],
  sp_dir[names(sp_dir) != "JORFLO"],
  linear_stack_predict, 
    sim = stack.iter,
    sum.fun = "median"
  )

### Model R-squared  ###
r2_list <- lapply(sp,function(x){
  
  # Subset data
  data <- fish_df %>% 
    left_join(
      sample_bridge,
      by = join_by(wateryear, region, site, species)
    ) %>% 
    filter(species == x) %>% 
    select(sample_id,length,age)
  
  # Candidate model R2
  ind_r2 <- len_R2(
    stack.df = sp_stack_wt[[x]],
    mod.dir = sp_dir[[x]],
    data = data,
    stack = F,
    sim = stack.iter,
    sum.fun = "median"
  )
  
  # Stack model R2
  stack_r2 <- len_R2(
    stack.df = sp_stack_wt[[x]],
    mod.dir = sp_dir[[x]],
    data = data,
    stack = T,
    sim = stack.iter,
    sum.fun = "median"
  )
  
  # bind into one data.frame
  rbind(ind_r2,stack_r2)
  
})
names(r2_list) <- sp


# Format predictions into data frames  -----------------------------------------

# Add species names to data frames
ind_mu_curve_list <- purrr::map2(
  ind_mu_curve_list,names(ind_mu_curve_list),
  function(x,y) x %>% mutate(species =y)
)
site_curve_list <- purrr::map2(
  site_curve_list,names(site_curve_list), 
  function(x,y) x %>% mutate(species =y)
  )
mu_curve_list <- purrr::map2(
  mu_curve_list,names(mu_curve_list), 
  function(x,y) x %>% mutate(species =y)
  )
ind_mean_growth_list <- purrr::map2(
  ind_mean_growth_list,names(ind_mean_growth_list), 
  function(x,y) x %>% mutate(species =y)
)
mean_growth_list <- purrr::map2(
  mean_growth_list,names(mean_growth_list), 
  function(x,y) x %>% mutate(species =y)
)
param_list <- purrr::map2(
  param_list,names(param_list), 
  function(x,y) x %>% mutate(species =y)
)
pred_list <- purrr::map2(
  pred_list,names(pred_list),
  function(x,y) x %>% mutate(species =y)
)
r2_list <- purrr::map2(
  r2_list,names(r2_list), 
  function(x,y) x %>% mutate(species =y)
  )

# Create one data.frame for all species
ind_mu_curve_df <- bind_rows(ind_mu_curve_list)
site_curve_df <- bind_rows(site_curve_list)
mu_curve_df <- bind_rows(mu_curve_list)
mean_growth_df <- bind_rows(mean_growth_list)
ind_mean_growth_df <- bind_rows(ind_mean_growth_list)
param_df <- bind_rows(param_list)
pred_df <- bind_rows(pred_list)
r2_df <- bind_rows(r2_list)


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
site_curve_bridged <-site_curve_df %>% 
  left_join(
    sample_bridge,
    by = join_by(species,sample_id)
    )


# Export  ----------------------------------------------------------------------

# For summary stats, export loo, stacking weights, growth rates, r2, 
# and parameters
saveRDS(
  sp_loo_compare,
  file.path(export_dir,paste0("loo_out_",Sys.Date(),".rds"))
  )
saveRDS(
  sp_stack_wt,
  file.path(export_dir,paste0("stack_wt_out_",Sys.Date(),".rds"))
  )
saveRDS(
  mean_growth_df, 
  file.path(
    export_dir,
    paste0("stacked_mean_growth_predictions_",Sys.Date(),".rds")
    )
  )
saveRDS(
  ind_mean_growth_df, 
  file.path(
    export_dir,
    paste0("ind_mean_growth_predictions_",Sys.Date(),".rds")
    )
  )
saveRDS(
  r2_df, 
  file.path(export_dir,paste0("model_r2_",Sys.Date(),".rds"))
  )
saveRDS(
  param_df, 
  file.path(export_dir,paste0("stacked_mu_parameters_",Sys.Date(),".rds"))
  )

# For plots, export length and growth-at-age predictions and bridge tables
saveRDS(
  site_curve_bridged, 
  file.path(export_dir,paste0("stacked_curves_",Sys.Date(),".rds"))
  )
saveRDS(
  mu_curve_df, 
  file.path(export_dir,paste0("stacked_mu_curves_",Sys.Date(),".rds"))
  )
saveRDS(
  ind_mu_curve_df, 
  file.path(export_dir,paste0("ind_mu_curves_",Sys.Date(),".rds"))
  )
saveRDS(
  pred_bridged, 
  file.path(
    export_dir,
    paste0("stacked_growth_predictions_",Sys.Date(),".rds")
    )
  )

