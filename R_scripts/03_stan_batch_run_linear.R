#-------------------------------------------------------------------------------
#
#   Stan Growth Curve Batch Processes and Diagnosis
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis
# CREATED: Feb 2, 2026

# DESCRIPTION: Prepares age-length and predictor data for rstan functions. Runs
# all stan models for each species, exporting results to species specific 
# directory, and return sampling diagnostics.


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load packages
library(dplyr)
library(readxl)
library(rstan)

# Directories
fun_dir <-"functions"
len_dir <- "~/Documents/Work/Everglades post-doc/Data analysis/Data cleaning/cleaned_data"
input_dir <- "input_data"
plot_dir <- "stan_outputs/plotting_info"
out_dir <- "stan_outputs/model_out_linear"

# Load in custom functions
source(file.path(fun_dir,"stan_loo_batch_functions.R"))

# Data (Make sure up-to date version!)
age_df <- readRDS(file.path(input_dir,"fsage_cleaned_2026-02-26.rds"))
len_df <- readRDS(file.path(len_dir,"fslen_cleaned_2026-02-25.rds"))
pred_df <-readRDS(file.path(input_dir,"fsgrw_predictors_2026-06-17.rds"))


# Create inputs for batch model runs  ------------------------------------------

# Create Stan data lists and id bridge tables for each species.
# Data is provide to Stan in named list for all variables. For random effects,
# groupings must be numeric. This function will create sampling event ids based
# on wateryear, region, and site. To link the new sampling id back to original
# site information, a bridge data.frame is also created
sp <- unique(age_df$species)

prep_list <-lapply(
  sp,
  stan_data_prep,
  age.df = age_df,
  len.df = len_df,
  pred.df=pred_df,
  fixed.effect = "linear",
  predictors = c("PC1","PC2","PC3"),
  scale = T,
  linear.predictions = T,
  pred.len = 100)
names(prep_list) <- sp

# Split data prep list into Stan data and plotting data (e.g. sample id, average
# length, and prediction labels) 
prep_list_t <- purrr::list_transpose(prep_list)
id_bridge <- bind_rows(prep_list_t$id_bridge)  # links sample_id to site and year
pred_lables <- prep_list_t$prediction_labels  # PC values used for predicted growth rates 
mean_lengths <- prep_list_t$mean_length  # mean length used to estimate inst. growth
input_list <- prep_list_t$stan_data   # data for stan analysis

# Export plotting data
saveRDS(id_bridge,file.path(plot_dir,paste0("fsgwh_sampleid_bridge_",Sys.Date(),".rds")))
saveRDS(pred_lables,file.path(plot_dir,paste0("fsgwh_pred_labels_",Sys.Date(),".rds")))
saveRDS(mean_lengths,file.path(plot_dir,paste0("fsgwh_mean_lengths_",Sys.Date(),".rds")))


# LUCGOO model runs  -----------------------------------------------------------

# Subset by species
luc <- "LUCGOO"
luc_input <- input_list[[luc]]

# Species specific nu for student t's distribution
luc_input$NU <- 4

# Run model, export results, and return diagnostics
luc_diag <- stan_diag_batch(
  fixed.effects = "linear",
  data = luc_input,
  sp = luc,
  export.dir = out_dir,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
  )


# POELAT model runs  -----------------------------------------------------------

# Subset by species
poe <- "POELAT"
poe_input <- input_list[[poe]]

# Species specific nu for student t's distribution
poe_input$NU <- 3

# Run model, export results, and return diagnostics
poe_diag <- stan_diag_batch(
  fixed.effects = "linear",
  data = poe_input,
  sp = poe,
  export.dir = out_dir,
  iter = 7000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97), 
  cores = 4
)
poe_diag_t <- purrr::transpose(poe_diag)


# HETFOR model runs  -----------------------------------------------------------

# Subset by species
het <- "HETFOR"
het_input <- input_list[[het]]

# Species specific nu for student t's distribution
het_input$NU <- 4

# Run model, export results, and return diagnostics
het_diag <- stan_diag_batch(
  fixed.effects = "linear",
  data = het_input,
  sp = het,
  export.dir = out_dir,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97), 
  mc.cores = 1
)
het_diag_t <- purrr::transpose(het_diag)


# GAMHOL model runs  -----------------------------------------------------------

# Subset by species
gam <- "GAMHOL"
gam_input <- input_list[[gam]]

# Species specific nu for student t's distribution
gam_input$NU <- 4

# Run model, export results, and return diagnostics
gam_diag <- stan_diag_batch(
  fixed.effects = "linear",
  data = gam_input,
  sp = gam,
  export.dir = out_dir,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)
gam_diag_t <- purrr::transpose(gam_diag)


# FUNCHR model runs  -----------------------------------------------------------

# Subset by species
fun <- "FUNCHR"
fun_input <- input_list[[fun]]

# Species specific nu for student t's distribution
fun_input$NU <- 3

# Run model, export results, and return diagnostics
fun_diag <- stan_diag_batch(
  fixed.effects = "linear",
  data = fun_input,
  sp = fun,
  export.dir = out_dir,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97),  
  cores = 4
)
fun_diag_t <- purrr::transpose(fun_diag)


# JORFLO model runs  -----------------------------------------------------------

# Subset by species
jor <- "JORFLO"
jor_input <- input_list[[jor]]

# Species specific nu for student t's distribution
jor_input$NU <- 4

# JORFLO has insufficient sample size for second level effects, run
# random effect only models
# Run model, export results, and return diagnostics
jor_diag <- stan_diag_batch(
  fixed.effects = NULL,
  data = jor_input,
  sp = jor,
  export.dir = out_dir,
  iter = 5000,
  warmup = 1000,
  chains =4,
  control = list(adapt_delta = .97), 
  cores = 4
)
jor_diag_t <- purrr::transpose(jor_diag)

