#-------------------------------------------------------------------------------
#
#  Simulation-based Stan Model testing
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: 02-27-2026

# DESCRIPTION: Simulates length-at-age and predictor variables based on different
# growth forms and effect structures to asses if  hierarchical growth models
# return the correct values.


# Housekeeping  ----------------------------------------------------------------
rm(list=ls())

# Packages
library(rstan)
library(loo)

# Directories
fun_dir <- "model_testing"

# Custom functions
source(file.path(fun_dir,"mod_test_functions.R"))


# Simulation inputs  -----------------------------------------------------------

# Parameters to monitor
params <- c("sigma_length",
            "mu_Linf","mu_g","mu_t",
            "cat_Linf","cat_g","cat_t",
            "tau",
            "beta_Linf","beta_g","beta_t")


# Data strucuture - all input
input_list <- list(
  n_sites = 18,
  n_ages = 25,
  max_age = 200,
  sigma_length = 1
)

# Random effect only inputs
vb_input <- list(
  mu_Linf =32,
  mu_g = 0.03,
  mu_t = -9,
  tau = c(0.02,0.03,2),
  cor.Linf.g =-.49,
  cor.Linf.t = .05,
  cor.g.t = 0.2
)
gz_input <- list(
  mu_Linf =32,
  mu_g = 0.009,
  mu_t = 41,
  tau = c(0.02,0.05,0.02),
  cor.Linf.g =-.49,
  cor.Linf.t = .05,
  cor.g.t = 0.2
)
lg_input <- list(
  mu_Linf =32,
  mu_g = 0.03,
  mu_t = 41,
  tau = c(0.02,0.05,0.02),
  cor.Linf.g =-.49,
  cor.Linf.t = .05,
  cor.g.t = 0.2
)

# linear inputs
linear_input <- list(
  beta_Linf =c(-.09,0,.1),
  beta_g =c(.09,-.01,.08),
  beta_t =c(0,.03,-.08)
)

# Effect inputs
vb_cat_input <- list(
  cat_Linf =c(28,32,36),
  cat_g =c(.03,.01,.03),
  cat_t =c(-9,-2,-11),
  tau = c(0.02,0.03,2),
  cor.Linf.g =-.49,
  cor.Linf.t = .05,
  cor.g.t = 0.2
)
gz_cat_input <- list(
  cat_Linf =c(28,32,36),
  cat_g =c(.02,.019,.027),
  cat_t =c(45,47,40),
  tau = c(0.02,0.05,.02),
  cor.Linf.g =-.49,
  cor.Linf.t = .05,
  cor.g.t = 0.2
)
lg_cat_input <- list(
  cat_Linf =c(28,32,36),
  cat_g =c(.029,.033,.03),
  cat_t =c(45,47,40),
  tau = c(0.02,0.05,.02),
  cor.Linf.g =-.49,
  cor.Linf.t = .05,
  cor.g.t = 0.2
)


# Simulate data  ---------------------------------------------------------------

mod_comb <-expand.grid(c("vb","gz","lg"),c("none","linear","categorical"))
mod_comb <-expand.grid(c("vb","gz","lg"),c("categorical"))

sim_list <- purrr::map2(
  mod_comb$Var1,
  mod_comb$Var2,
  function(mod.form,fixed.effects){
  
    # Prepare effect structure-specific inputs
    if(fixed.effects == "none") {
      if(mod.form == "vb") input_list <- c(input_list,vb_input)
      if(mod.form == "gz") input_list <- c(input_list,gz_input)
      if(mod.form == "lg") input_list <- c(input_list,lg_input)
    }
    
    if(fixed.effects == "linear") {
      if(mod.form == "vb") input_list <- c(input_list,vb_input,linear_input)
      if(mod.form == "gz") input_list <- c(input_list,gz_input,linear_input)
      if(mod.form == "lg") input_list <- c(input_list,lg_input,linear_input)
    }
    
    if(fixed.effects == "categorical"){
      if(mod.form == "vb") input_list <- c(input_list,vb_cat_input)
      if(mod.form == "gz") input_list <- c(input_list,gz_cat_input)
      if(mod.form == "lg") input_list <- c(input_list,lg_cat_input)
    } 
    
    # Create param list and actual values
    param_select <- params[params %in% names(input_list)]
    actual_values <- unlist(input_list[param_select])
  
    # simulate data for specified model
    sim <-replicate(
      10,
      ageLengthSim(
        sim.input = input_list,
        nu=3,
        mod.form = mod.form,
        fixed.effects = fixed.effects,
        ),
      simplify = F
      )
    
    # Prepare output
    out <- list(sim,mod.form,fixed.effects,param_select,actual_values)
    names(out) <-c(
      "sim.list",
      "mod.form",
      "fixed.effects",
      "params",
      "actual.values"
      # equal.cat =T
      )
    out
    }
  )
names(sim_list) <- paste(mod_comb$Var1,mod_comb$Var2,sep = "_")
sim_list_t <- purrr::transpose(sim_list)


# Model runs  ------------------------------------------------------------------

# Fit each type of model to simulated data sets
sim_out <- Map(
  function(sim.list,mod.form,fixed.effect,params){
    
    # Run model with each simulate data
    parallel::mclapply(
      sim.list,
      model_sim_test,
        params = params,
        mod.form = mod.form,
        fixed.effects = fixed.effect,
        chains = 4,
        warmup=1000,
        iter=3000,
        cores=1,
        refresh = 0,
        control = list(adapt_delta = 0.97),
      mc.cores = 10
    )
  },
  sim_list_t$sim.list,
  sim_list_t$mod.form,
  sim_list_t$fixed.effects,
  sim_list_t$params
)


# Result summaries  ------------------------------------------------------------

purrr::map2(sim_out, sim_list_t$actual.values, function(sim.out,actual.values) {
  
  # Prepare list
  out_list <- purrr::transpose(sim.out)
  
  # Summarize mean results 
  sim_out_df <- do.call(rbind,out_list[[1]])
  sim_out_med <- apply(sim_out_df,2,quantile,.5,na.rm=T)
  sim_out_lwr <- apply(sim_out_df,2,quantile,.025,na.rm=T)
  sim_out_upr <- apply(sim_out_df,2,quantile,.975,na.rm=T)
  
  # Summarize overlap
  sim_overlap_df <- do.call(rbind,out_list[[2]])
  sim_overlap_sum <- colSums(sim_overlap_df)/nrow(sim_overlap_df)
  
  # Compare results  
  rbind(actual.values,sim_out_med,sim_out_lwr,sim_out_upr,sim_overlap_sum)
  
})

