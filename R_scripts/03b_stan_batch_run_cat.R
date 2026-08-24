#-------------------------------------------------------------------------------
#
#   Stan Growth Curve Batch Processes and Diagnosis
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis
# CREATED: Feb 2, 2026

# DESCRIPTION: Runs all stan categorical fixed effect models for each species 
# (except JORFLO), exporting results to species specific directory, and return 
# sampling diagnostics.


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load packages
library(dplyr)
library(readxl)
library(rstan)

# Directories
# fun_dir <-"functions"
len_dir <- paste0(
  "~/Documents/Work/Everglades post-doc/",
  "Data analysis/Data cleaning/cleaned_data"
  )
input_dir <- "input_data"
out_dir <- "outputs/stan_outputs"
fig_dir <-"figures"

# Load in custom functions
devtools::load_all("~/Documents/work/R packages/growthstack")

# Data (Make sure up-to date version!)
age_df <- readRDS(file.path(input_dir,"fsage_cleaned_2026-06-18.rds"))
len_df <- readRDS(file.path(len_dir,"fslen_cleaned_2026-02-25.rds"))
pred_df <-readRDS(file.path(input_dir,"fsgrw_predictors_2026-08-21.rds"))

# Combine age and predictor data.frames
input_df <- age_df %>% 
  left_join(pred_df)

# LUCGOO model runs  -----------------------------------------------------------

# Fit models
luc_out <- fit_growth(
  mod.form = c("vb","gz","lg"),
  nu=4,
  fixed.effect = "categorical",
  sample.groups = c("wateryear","region","site"),
  predictors = c("PC1","PC2","PC3"),  
  scale = T,
  linear.predictions = T,  
  pred.len = 100,
  sp="LUCGOO",   
  age.df = input_df, 
  len.df = len_df,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)

# Check for sampling and convergence issues
lapply(luc_out$model_out, stan_diag)

# export
luc_dir <- file.path(out_dir,"LUCGOO")
dir.create(luc_dir)
lapply(1:length(luc_out$model_out),function(x){
  mod <- luc_out$model_out[[x]]
  name <- names(luc_out$model_out)[x]
  file_name <- paste0(
    name,
    "_LUCGOO_",
    Sys.Date(),
    ".rds"
  )
  saveRDS(mod,file.path(luc_dir,file_name))
})


# POELAT model runs  -----------------------------------------------------------

# Fit models
poe_out <- fit_growth(
  mod.form = c("vb","gz","lg"),
  nu=3,
  fixed.effect = "categorical",
  sample.groups = c("wateryear","region","site"),
  predictors = c("PC1","PC2","PC3"),  
  scale = T,
  linear.predictions = T,  
  pred.len = 100,
  sp="POELAT",   
  age.df = input_df, 
  len.df = len_df,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)

# Check for sampling and convergence issues
lapply(poe_out$model_out, stan_diag)

# export
poe_dir <- file.path(out_dir,"POELAT")
dir.create(poe_dir)
lapply(1:length(poe_out$model_out),function(x){
  mod <- poe_out$model_out[[x]]
  name <- names(poe_out$model_out)[x]
  file_name <- paste0(
    name,
    "_POELAT_",
    Sys.Date(),
    ".rds"
  )
  saveRDS(mod,file.path(poe_dir,file_name))
})


# HETFOR model runs  -----------------------------------------------------------

# Fit models
het_out <- fit_growth(
  mod.form = c("vb","gz","lg"),
  nu=4,
  fixed.effect = "categorical",
  sample.groups = c("wateryear","region","site"),
  predictors = c("PC1","PC2","PC3"),  
  scale = T,
  linear.predictions = T,  
  pred.len = 100,
  sp="HETFOR",   
  age.df = input_df, 
  len.df = len_df,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)

# Check for sampling and convergence issues
lapply(het_out$model_out, stan_diag)

# export
het_dir <- file.path(out_dir,"HETFOR")
dir.create(het_dir)
lapply(1:length(het_out$model_out),function(x){
  mod <- het_out$model_out[[x]]
  name <- names(het_out$model_out)[x]
  file_name <- paste0(
    name,
    "_HETFOR_",
    Sys.Date(),
    ".rds"
  )
  saveRDS(mod,file.path(het_dir,file_name))
})


# GAMHOL model runs  -----------------------------------------------------------

# Fit models
gam_out <- fit_growth(
  mod.form = c("vb","gz","lg"),
  nu=4,
  fixed.effect = "categorical",
  sample.groups = c("wateryear","region","site"),
  predictors = c("PC1","PC2","PC3"),  
  scale = T,
  linear.predictions = T,  
  pred.len = 100,
  sp="GAMHOL",   
  age.df = input_df, 
  len.df = len_df,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)

# Check for sampling and convergence issues
lapply(gam_out$model_out, stan_diag)

# export
gam_dir <- file.path(out_dir,"GAMHOL")
dir.create(gam_dir)
lapply(1:length(gam_out$model_out),function(x){
  mod <- gam_out$model_out[[x]]
  name <- names(gam_out$model_out)[x]
  file_name <- paste0(
    name,
    "_GAMHOL_",
    Sys.Date(),
    ".rds"
  )
  saveRDS(mod,file.path(gam_dir,file_name))
})


# FUNCHR model runs  -----------------------------------------------------------

# Fit models
fun_out <- fit_growth(
  mod.form = c("vb","gz","lg"),
  nu=3,
  fixed.effect = "categorical",
  sample.groups = c("wateryear","region","site"),
  predictors = c("PC1","PC2","PC3"),  
  scale = T,
  linear.predictions = T,  
  pred.len = 100,
  sp="FUNCHR",   
  age.df = input_df, 
  len.df = len_df,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)

# Check for sampling and convergence issues
lapply(fun_out$model_out, stan_diag)

# export
fun_dir <- file.path(out_dir,"FUNCHR")
dir.create(fun_dir)
lapply(1:length(fun_out$model_out),function(x){
  mod <- fun_out$model_out[[x]]
  name <- names(fun_out$model_out)[x]
  file_name <- paste0(
    name,
    "_FUNCHR_",
    Sys.Date(),
    ".rds"
  )
  saveRDS(mod,file.path(fun_dir,file_name))
})


# Category labels  -------------------------------------------------------------

out_list <- list(
  LUCGOO = luc_out,
  POELAT = poe_out,
  HETFOR = het_out,
  GAMHOL = gam_out,
  FUNCHR = fun_out
)


# Split data prep various plotting data (e.g. sample id, average
# length, and prediction labels) 
out_list_t <- purrr::list_transpose(out_list)
id_bridge <- bind_rows(out_list_t$id_bridge)  # links sample_id to site and year
cat_lables <- bind_rows(out_list_t$category_labels)  # hydroperiod classes


# Export plotting data]
saveRDS(
  cat_lables,
  file.path(
    fig_dir,
    "_labels",
    paste0("fsgwh_cat_labels",Sys.Date(),".rds")
    )
  )

