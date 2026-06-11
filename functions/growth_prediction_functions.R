# grouped prediction df function  ----------------------------------------------

# DESCRIPTION: Prepares input dataframe for length2age_predict and  
# length2growth_predict functions for population, sampling event, and hydroperiod
# level predictions. Also produces bridge function to link combined grouping ids
# back to original group level ids.

# INPUT:
# ! group.size =vector containing the number of groups in the population (1),
#                sampling event, and hydroperiod groups
# ! sp = character to indicate which species is being modelled. Used in bridge
#        dataframe
# ! min.pred = minimum prediction input value
# ! max.pred = maximum prediction input value

# OUTPUT: a named list containing the prediction function input dataframe and
# a bridge dataframe to link combined grouping ids to sampling and hydroperiod
# ids. Prediction dataframe has a column for the combined grouping id (group_id) 
# and the prediction input value (pred; length or age)

# REQUIRES: dplyr (all), tidyr,

grouping_predDF <- function(group.size,sp,min.pred,max.pred,group_vec = c("mu_","hydr_","")) {
  
  # # Three groupings used in model
  # group_vec <- c("mu_","hydr_","")
  
  # Create a unique id for each grouping among all groups
  id_df <- data.frame(group = unlist(purrr::map2(group_vec,group.size,rep)),
                      old_id = unlist(purrr::map(group.size,seq,from=1)),
                      group_id=1:sum(group.size))
  
  # Create prediction dataframe with a user specified range of input (age or length)
  # data for each new grouping. This can be feed to one of the prediction
  # functions
  pred_df <- id_df %>% 
    tidyr::crossing(input = min.pred:max.pred) %>% 
    select(group_id,input)
  
  # Create a dataframe that links to new group ids to hydroperiod or sample
  # event ids. This wide format table can be merged directly into the 
  # sample_id_bridge dataframe
  id_bridge_df<-  id_df %>% 
    mutate(
      species = sp,  
      group = case_when(
        group == "" ~ "sample_id",
        group == "hydr_" ~ "hydro_id",
        group == "mu_" ~ "mu",
        T ~ NA
      )) %>% 
    tidyr::pivot_wider(names_from = group,
                       values_from = old_id)
  
  # return a named list
  out <- list(pred_df,id_bridge_df)
  names(out) <- c("prediction","id_bridge")
  out
}


# group predict stacking  ------------------------------------------------------

# DESCRIPTION: Create model stacked, group specific growth curve predictions or 
# growth parameters estimates. Function extract parameter posterior distributions 
# based on model names provided with stacking weights. Wrapper for param_extract, 
# prediction_sampler, and parameter_sampler

# INPUT:
# ! stack.df = dataframe containing model file names and stacking wts. Must have
#              columns "model" and "stack_wt
# ! mod.dir  = pathway for Stan model output files
# ! group.id = Character("hydr_","mu_", "","all") indicating the grouping level  
#               of model parameters to extract. If all is selected, population,
#               site, and hydroperiod level parameters will be extracted
# ! group.size = number of groups in grouping level# ! sim = total number of postieror draws for stacking
# ! type = Character ("parameter" or "prediction"), model stack growth predictions
#          or paramters
# ! ... = additional arguments for the prediction_sampler and parameter_sampler
#         functions

# OUTPUT: 3D array containing either group specific model parameters, or predictions
# with each slcie being an random psoterior draw from one of the models in stack.

# REQUIRES: abind

growth_stackR <- function(stack.df, mod.dir, group.id, group.size, sim,type,sum.fun, ...){
  
  # Prediction or parameters?
  fun_name = paste0(type,"_sampler")
  fun = get(fun_name)
  
  # Retain models with 
  stack <- stack.df %>% 
    mutate(n_sim = round(sim*stack_wt)) %>% 
    filter(n_sim >0)
  
  # Load it model parameters
  mods<- stack$model
  
  # Select desired groupings
  if (group.id == "all") {
    stopifnot("Please provide group sizes for all three groupings" = 
                length(group.size) == 3)
    mod_list <- lapply(mods,param_extract_all,
                       group.size = group.size,
                       mod.dir = mod.dir)
  } else {
    mod_list <- lapply(mods,param_extract,mod.dir,group.id,group.size)
  }
  
  # prepare prediction inputs
  sim_list = stack$n_sim
  g_mod_list <- substr(mods,1,2)
  
  # Run prediction function
  pred_list <-Map(function(x,y,z) fun(model.out = x,
                                      g.mod = y, 
                                      n.sim = z, 
                                      ...=...),
                  mod_list, g_mod_list, sim_list)
  
  # Bind results into 3d array
  out <- abind::abind(pred_list, along = 3)
  
  # Set grouping for summary function
  args <- list(...)
  if("input.var" %in% names(args)){
    group_var <- c("group_id",args$input.var)
  } else {
    group_var <- c("group_id")
  }
  
  # Summarize into dataframe
  boot_summary(
    out,
    sum.fun=sum.fun,
    group.var=group_var
  )
  
}

# length predictions of all models  --------------------------------------------

# DESCRIPTION: Estiame R2 of length predictions of all models in stacking output. 
# Function extract parameter posterior distributions and predicts teh length at
# age for each grouping. R2 is then estimated for each model Wrapper for , 
# param_extract and prediction_sampler.

# INPUT:
# ! stack.df = dataframe containing model file names and stacking wts. Must have
#              columns "model" and "stack_wt
# ! mod.dir  = pathway for Stan model output files
# ! group.id = Character("hydr_","mu_", "","all") indicating the grouping level  
#               of model parameters to extract. If all is selected, population,
#               site, and hydroperiod level parameters will be extracted
# ! group.size = number of groups in grouping level# ! sim = total number of 
#                postieror draws for stacking
# ! n.sims = number of posteiror draws to estimate mean with
# ! ... = additional arguments for the prediction_sampler functions

# OUTPUT: dataframe containing r2 and adjusted r2 for each model

# REQUIRES: purrr

curve_predictR <- function(stack.df, mod.dir, group.id, group.size,n.sim, type,sum.fun="mean", ...){
  
  # Prediction or parameters?
  fun_name = paste0(type,"_sampler")
  fun = get(fun_name)
  
  # Load it model parameters
  mods<- stack.df$model
  
  # Extract posterior distributions
  mod_list <- lapply(mods,param_extract,mod.dir,group.id,group.size)
  
  # prepare prediction inputs
  g_mod_list <- substr(mods,1,2)
  
  # Run prediction function
  pred_list <-Map(function(x,y) fun(model.out = x,
                                    g.mod = y, 
                                    n.sim = n.sim,
                                    ...),
                  mod_list, g_mod_list)
  names(pred_list) <- mods
  
  # Set grouping for summary function
  args <- list(...)
  if("input.var" %in% names(args)){
    group_var <- c("group_id",args$input.var)
  } else {
    group_var <- c("group_id")
  }
  
  # Extract means
  pred_summary <- lapply(pred_list, boot_summary,sum.fun=sum.fun, group.var=group_var)
  
  # Add model name
  pred_summary <- purrr::map2(pred_summary,names(pred_summary),
                              function(x,y) x %>% mutate(mod =y))
  
  # Combine into one data frame
  bind_rows(pred_summary)
  
}


# Estimate R2 of length predictions of all models  ------------------------------

# DESCRIPTION: Estiame R2 of length predictions of all models in stacking output. 
# Function extract parameter posterior distributions and predicts teh length at
# age for each grouping. R2 is then estimated for each model Wrapper for , 
# param_extract and prediction_sampler.

# INPUT:
# ! stack.df = dataframe containing model file names and stacking wts. Must have
#              columns "model" and "stack_wt
# ! mod.dir  = pathway for Stan model output files
# ! group.id = Character("hydr_","mu_", "","all") indicating the grouping level  
#               of model parameters to extract. If all is selected, population,
#               site, and hydroperiod level parameters will be extracted
# ! group.size = number of groups in grouping level# ! sim = total number of postieror draws for stacking
# ! n.sims = number of posteiror draws to estimate mean with
# ! ... = additional arguments for the prediction_sampler functions

# OUTPUT: dataframe containing r2 and adjusted r2 for each model

# REQUIRES: purrr


len_R2 <- function(stack.df, mod.dir, group.id, group.size, n.sim,sp,sum.fun="mean", ...){
  
  # Load it model parameters
  mods<- stack.df$model
  
  # Extract posterior distributions
  mod_list <- lapply(mods,param_extract,mod.dir,group.id,group.size)
  
  # prepare prediction inputs
  g_mod_list <- substr(mods,1,2)
  
  # Run prediction function
  pred_list <-Map(function(x,y) prediction_sampler(model.out = x,
                                                   g.mod = y, 
                                                   n.sim = n.sim,
                                                   input.var = "age",
                                                   output.var = "length",
                                                   ...),
                  mod_list, g_mod_list)
  names(pred_list) <- mods
  
  # Extract means
  pred_summary <- lapply(pred_list, boot_summary,sum.fun = sum.fun, group.var=c("group_id","age"))
  
  # Add model name
  pred_summary <- purrr::map2(pred_summary,names(pred_summary),
                              function(x,y) x %>% mutate(mod =y))
  
  # Combine into one data frame
  pred_df <- bind_rows(pred_summary)
  
  # Connect group id to site info
  bridge_df <- curve_id_bridge %>% 
    filter (species == sp) %>% 
    left_join(sample_bridge, by = join_by(species,sample_id)) %>% 
    right_join(pred_df,by = join_by(group_id))%>%
    select(-group_id,-sample_id)
  
  # Link to actual values
  actual_df <- fish_df %>% 
    filter(species == sp) %>% 
    rename(age = ring_count) %>% 
    left_join(bridge_df, by = join_by(wateryear,region,site,age,species))
  
  # Calculate r2 values of predictions for each model
  r2_list <- lapply(mods, function (x) {
    
    # subset data for select mode
    df <- actual_df %>% 
      filter(mod == x)
    
    # Estimate r2
    out <- lm(length_pred_median ~ length,df)
    data.frame(model = x,r2 = summary(out)$r.squared, adj_r2 = summary(out)$adj.r.squared)
  }) 
  bind_rows(r2_list)
}


# linear predict stacking  -----------------------------------------------------
  
# DESCRIPTION: Create model stacked, group specific parameter and inst. growth
# predictions. Function extract prediction posterior distributions based on model 
# names provided with stacking weights. Wrapper for predict_extract, 
# linear_prediction_sampler.

# INPUT:
# ! stack.df = dataframe containing model file names and stacking wts. Must have
#              columns "model" and "stack_wt
# ! mod.dir  = pathway for Stan model output files

# OUTPUT: 3D array containing  predictions with each slcie being an random 
# posterior draw from one of the models in stack.

# REQUIRES: abind

linear_pred_stackR <- function(stack.df, mod.dir, sim,sum.fun, ...){
  
  # Retain models with 
  stack <- stack.df %>% 
    mutate(n_sim = round(sim*stack_wt)) %>% 
    filter(n_sim >0)
  
  # Load it model parameters
  mods<- stack$model
  sim_list = stack$n_sim
  
  # Extract posterior distributions of all models with stacking weights
  mod_list <- lapply(mods,predict_extract,mod.dir)
  
  # Run prediction function
  pred_list <- purrr::map2(mod_list,sim_list, linear_predict_sampler)
  
  # Bind results into 3d array
  out <- abind::abind(pred_list, along = 3)
  
  # Summarize into dataframe
  boot_summary(
    out,
    sum.fun=sum.fun,
    group.var=c("pred_id")
  )
  
}  
