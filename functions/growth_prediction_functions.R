#-------------------------------------------------------------------------------
#
# Growth curve and parameter prediction functions
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: Feb 5, 2026

# DESCRIPTION: Functions used to create both candidate and model-stacked growth
# curve predictions using Stan growth model outputs. Functions also can model-
# stack growth parameters across models. Script includes helper functions used
# in main functions


# grouping_predDF  -------------------------------------------------------------

#' Grouped prediction data.frame preparation.
#' 
#' @description Prepares input data.frame for length2age_predict and  
#' length2growth_predict functions for population, sampling event, and hydroperiod
#' level predictions. Also produces bridge function to link combined grouping ids
#' back to original group level ids.
#'
#' @param group.size Vector containing the number of groups in the population 
#' (1),sampling event, and hydroperiod groups.
#' @param sp Character to indicate which species is being modeled. Used in 
#' bridgedata.frame
#' @param min.pred Minimum prediction input value
#' @param max.pred Maximum prediction input value
#' 
#' 
#' @returns a named list containing the prediction function input data.frame and 
#' a bridge data.frame to link combined grouping ids to sampling and hydroperiod 
#' ids. Prediction data.frame has a column for the combined grouping id 
#' (group_id) and the prediction input value (pred; length or age).
#' @export

# REQUIRES: dplyr (all), tidyr,

grouping_predDF <- function(group.size,sp,min.pred,max.pred,group_vec = c("mu_","hydr_","")) {

  # Create a unique id for each grouping among all groups
  id_df <- data.frame(group = unlist(purrr::map2(group_vec,group.size,rep)),
                      old_id = unlist(purrr::map(group.size,seq,from=1)),
                      group_id=1:sum(group.size))
  
  # Create prediction data.frame with a user specified range of input (age or length)
  # data for each new grouping. This can be feed to one of the prediction
  # functions
  pred_df <- id_df %>% 
    tidyr::crossing(input = min.pred:max.pred) %>% 
    select(group_id,input)
  
  # Create a data.frame that links to new group ids to hydroperiod or sample
  # event ids. This wide format table can be merged directly into the 
  # sample_id_bridge data.frame
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


# growth_stackR  ---------------------------------------------------------------

#' Model stack predictions and growth parameters
#' 
#' @description Create model stacked, group specific growth curve predictions 
#' or growth parameters estimates with credible intervals. 
#' 
#' @param stack.df Data.frame containing model file names and stacking wts. 
#' Must have columns "model" and "stack_wt
#' @param mod.dir  File path for Stan model output files
#' @param group.id Character ("group","mu", "site") indicating the grouping 
#' level of model parameters to extract. 
#' @param sim Number of posterior draws for stacking
#' @param type Character ("parameter" or "prediction"), model stack growth 
#' predictions or parameters
#' @param type sum.fun Character ("mean" or "median) for type of summary 
#' statistic of posterior distribution. Default is mean 
#' @param ... = Additional arguments passed to .prediction_sampler or 
#' .parameter_sampler auxiliary functions. Required arguments for stacked growth 
#' curve predictions include:
#'   \describe{
#'     \item{input.df}{Data.frame with a column ("input) for prediction input  
#'      data (length or age) and an optional column ("group_id") for a grouping 
#'     variable}
#'     \item{input.var}{Character ("length" or "age") to indicate if prediction 
#'     is made using length or age data}
#'     \item{output.var}{Character ("length","age","growth","interval_growth") 
#'     to indicate what type of prediction is being made. If vector is provided, 
#'     multiple prediction columns will be created. Growth indicates 
#'     instantaneous growth, and interval growth is exponential growth during a 
#'     specified interval}
#'     \item{days}{Number of days to estimate interval growth. Only required 
#'     if output.var includes "interval_growth"}
#'     \item{wt.df}{wt.df = Data.frame containing length-weight parameters. 
#'     Only required if output.var includes "interval_growth"}
#'     \item{dry.wt}{Optional dry weight conversion factor. Only required 
#'     if output.var includes "interval_growth". Default = 1}
#'     \item{parallel}{T or F. Use multiple cores. Only should be used on 
#'     Linux and MacOS. Default is F.}
#'     \item{mc.cores}{Number of core for parallel processing if parallel = T
#'     Default is NULL}
#'   }
#'   
#'Required arguments for stacked parameter estimates include:
#'   \describe{
#'     \item{truncate.inf}{Truncate negative values of the inflection 
#'     parameter to zero (i.e., faster growth at birth). Default is F. }
#'   }
#'   
#' @details Function extract parameter posterior distributions based on model  
#' names provided with stacking weights. Combines posterior distributions for 
#' predictions based on model stacking weights. Summarizes posterior 
#' distributions for predictions with mean or median, and 95% credible intervals. 
#' Wrapper for .param_extract, .prediction_sampler, .parameter_sampler, and 
#' .boot_summary
#' 
#' @returns Data.frame with a row for every input (age or length) and grouping
#' combination. Contains columns for age or length, grouping index, and columns 
#' for prediction summary statistics (mean or median, lower, and upper 
#' credible interval) for each output variable.
#' 
#' @export

# REQUIRES: abind

growth_stackR <- function(stack.df, mod.dir, group.id, sim,type,sum.fun, ...){
  
  # Prediction or parameters?
  fun_name = paste0(".",type,"_sampler")
  fun = get(fun_name)
  
  # Retain models with 
  stack <- stack.df %>% 
    mutate(n_sim = round(sim*stack_wt)) %>% 
    filter(n_sim >0)
  
  # Load it model parameters
  mods<- stack$model
  
  # Select desired groupings
  mod_list <- lapply(mods,.param_extract,mod.dir,group.id)

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
  
  # Summarize into data.frame
  .boot_summary(
    out,
    sum.fun=sum.fun,
    group.var=group_var
  )
  
}


# curve_predictR  --------------------------------------------------------------

#' Candidate model predictions and growth parameters
#' 
#' @description Create group specific growth curve predictions 
#' or growth parameters estimates with credible intervals for each candidate
#' model used in model stacking. 
#' 
#' 
#' @param stack.df Data.frame containing model file names and stacking wts. 
#' Must have columns "model" and "stack_wt
#' @param mod.dir  File path for Stan model output files
#' @param group.id Character ("group","mu", "site") indicating the grouping 
#' level of model parameters to extract.
#' @param n.sim Number of posterior draws for stacking
#' @param type Character ("parameter" or "prediction"), model stack growth 
#' predictions or parameters
#' @param type sum.fun Character ("mean" or "median) for type of summary 
#' statistic of posterior distribution. Default is mean 
#' @param ... = Additional arguments passed to .prediction_sampler or 
#' .parameter_sampler auxiliary functions. Required arguments for stacked growth 
#' curve predictions include:
#'   \describe{
#'     \item{input.df}{Data.frame with a column ("input) for prediction input  
#'      data (length or age) and an optional column ("group_id") for a grouping 
#'     variable}
#'     \item{input.var}{Character ("length" or "age") to indicate if prediction 
#'     is made using length or age data}
#'     \item{output.var}{Character ("length","age","growth","interval_growth") 
#'     to indicate what type of prediction is being made. If vector is provided, 
#'     multiple prediction columns will be created. Growth indicates 
#'     instantaneous growth, and interval growth is exponential growth during a 
#'     specified interval}
#'     \item{days}{Number of datys to estimate interval growth. Only required 
#'     if output.var includes "interval_growth"}
#'     \item{wt.df}{wt.df = Data.frame containing length-weight parameters. 
#'     Only required if output.var includes "interval_growth"}
#'     \item{dry.wt}{Optional dry weight conversion factor. Only required 
#'     if output.var includes "interval_growth". Default = 1}
#'     \item{parallel}{T or F. Use multiple cores. Only should be used on 
#'     Linux and MacOS. Default is F.}
#'     \item{mc.cores}{Number of core for parallel processing if parallel = T
#'     Default is NULL}
#'   }
#'   
#'Required arguments for stacked parameter estimates include:
#'   \describe{
#'     \item{truncate.inf}{Truncate negative values of the inflection 
#'     parameter to zero (i.e., faster growth at birth). Default is F }
#'   }
#'   
#' @details Function extract parameter posterior distributions based on model  
#' names provided with stacking weights. Summarizes posterior distributions for 
#' predictions with mean or median, and 95% credible intervals for each 
#' candidate model. Wrapper for .param_extract, .prediction_sampler, 
#' .parameter_sampler, and .boot_summary
#' 
#' @returns Data.frame with a row for every input (age or length), grouping,
#' and model combination. Contains columns for age or length, model name, 
#' grouping index, and columns for prediction summary statistics (mean or 
#' median, lower-lwr, and upper-upr credible interval) for each output variable.
#' 
#' @export

# REQUIRES: purrr

curve_predictR <- function(stack.df, mod.dir, group.id,n.sim, type,sum.fun="mean", ...){
  
  # Prediction or parameters?
  fun_name = paste0(".",type,"_sampler")
  fun = get(fun_name)
  
  # Load it model parameters
  mods<- stack.df$model
  
  # Extract posterior distributions
  mod_list <- lapply(mods,.param_extract,mod.dir,group.id)
  
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
  pred_summary <- lapply(pred_list, .boot_summary,sum.fun=sum.fun, group.var=group_var)
  
  # Add model name
  pred_summary <- purrr::map2(pred_summary,names(pred_summary),
                              function(x,y) x %>% mutate(mod =y))
  
  # Combine into one data.frame
  bind_rows(pred_summary)
  
}


# len_R2  ----------------------------------------------------------------------

#' R-sqaured estimation for candidate models.
#' 
#' @description Estimates r-squared of length predictions of all candidate 
#' models in stacking output. 
#' 
#' @param stack.df Data.frame containing model file names and stacking wts. 
#' Must have columns "model" and "stack_wt
#' @param mod.dir File path for Stan model output files
#' @param group.id Character ("group","mu", "site") indicating the grouping 
#' level of model parameters to extract.
#' @param n.sims Number of posterior draws
#' @param sp Character for name of species of which to estimate r-squared
#' @param sum.fun Character ("mean" or "median) for type of summary 
#' statistic of posterior distribution. Default is mean 
#' @param ... = Additional arguments passed to .prediction_sampler or 
#' .parameter_sampler auxiliary functions. Required arguments for stacked growth 
#' curve predictions include:
#'   \describe{
#'     \item{input.df}{Data.frame with a column ("input) for prediction input  
#'      data (length or age) and an optional column ("group_id") for a grouping 
#'     variable}
#'     \item{parallel}{T or F. Use multiple cores. Only should be used on 
#'     Linux and MacOS. Default is F.}
#'     \item{mc.cores}{Number of core for parallel processing if parallel = T.
#'     Default is NULL.}
#'   }
#'
#' 
#' @details Function extract parameter posterior distributions for each model and 
#' predicts the length at age for each grouping. R-squared is then estimated 
#' for each model. Wrapper for .param_extract and .prediction_sampler.
#' 
#' @returns Data.frame containing r-squared and adjusted r-squared for each 
#' model
#' @export

# REQUIRES: purrr

len_R2 <- function(stack.df, mod.dir, group.id, n.sim,sp,sum.fun="mean", ...){
  
  # Load it model parameters
  mods<- stack.df$model
  
  # Extract posterior distributions
  mod_list <- lapply(mods,.param_extract,mod.dir,group.id)
  
  # prepare prediction inputs
  g_mod_list <- substr(mods,1,2)
  
  # Run prediction function
  pred_list <-Map(function(x,y) 
    .prediction_sampler(model.out = x,
       g.mod = y, 
       n.sim = n.sim,
       input.var = "age",
       output.var = "length",
       ...),
    mod_list, 
    g_mod_list)
  names(pred_list) <- mods
  
  # Extract means
  pred_summary <- lapply(pred_list, .boot_summary,sum.fun = sum.fun, group.var=c("group_id","age"))
  
  # Add model name
  pred_summary <- purrr::map2(pred_summary,names(pred_summary),
                              function(x,y) x %>% mutate(mod =y))
  
  # Combine into one data.frame
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


# linear_pred_stackR  ----------------------------------------------------------
  
#' Model-stacked instantaneous growth predictions
#'
#' @description Creates model-stacked predictions of instantaneous growth rate
#' across a range of predictor variables.
#' 
#' @param stack.df Data.frame containing model file names and stacking weights. 
#' Must have columns "model" and "stack_wt. Models must include predictions of
#' instantaneous growth (inst_growth) and asymptotic length (Linf) across a 
#' range of predictor variables.
#' @param mod.dir  File path for Stan model output files
#' @param sim Number of posterior draws for stacking
#' @param sum.fun Character ("mean" or "median) for type of summary 
#' statistic of posterior distribution. Default is mean 
#'
#' @details Function extract prediction posterior distributions of instantaneous 
#' growth rates across predictors based on model names provided with stacking 
#' weights. Posterior distributions are combined based on stacking weights and 
#' summarized with mean or median, and 95% credible intervals. Wrapper for 
#' .predict_extract, linear_prediction_sampler.
#' 
#' @returns Data.frame with rows for each prediction input (e.g. PC1 value) and
#' columns for prediction number (pred_id), and columns for prediction summary 
#' statistics (mean or median, lower-lwr, and upper-upr credible interval) for 
#' each output variable.
#' 
#' @export

# REQUIRES: abind

linear_pred_stackR <- function(stack.df, mod.dir, sim,sum.fun){
  
  # Retain models with 
  stack <- stack.df %>% 
    mutate(n_sim = round(sim*stack_wt)) %>% 
    filter(n_sim >0)
  
  # Load it model parameters
  mods<- stack$model
  sim_list = stack$n_sim
  
  # Extract posterior distributions of all models with stacking weights
  mod_list <- lapply(mods,.predict_extract,mod.dir)
  
  # Run prediction function
  pred_list <- purrr::map2(mod_list,sim_list, .linear_predict_sampler)
  
  # Bind results into 3d array
  out <- abind::abind(pred_list, along = 3)
  
  # Summarize into data.frame
  .boot_summary(
    out,
    sum.fun=sum.fun,
    group.var=c("pred_id")
  )
  
}  

# Model parameter extraction helpers -------------------------------------------

# Extracts posterior distribution of the asymptotic, scaling, and inflection
# parameters from a growth model stanfit object based on user selected
# groupings (grouping - groupt, sampling event - site, or population - mu).
# Returns a 3d array (num groups, 4, iterations). Each slice is a posterior
# draw of the asymptotic, scaling, and inflection parameter values for each
#  group

# REQUIRES: rstan, abind

.param_extract <-function(mod.out,mod.dir,group.id){
  
  # Each model has different name for same class of parameter, create list with
  # all the names for each class
  param.list <- list(
    Linf = c("Linf"),
    slope = c("g1","g2","g3"),
    inf = c("t0","ti")
  )
  
  # Load in model outputs
  mod <- readRDS(file.path(mod.dir,mod.out))
  
  # Load in posterior samples for all parameters
  sim.list <- rstan::extract(mod,permute =T)
  
  # How many posterior draws are there?
  n.iter <- dim(sim.list[[1]])[1]
  
  # Extract draws from group-specific group parametes
  out.list <-lapply(param.list, 
                    .single_param_extract,
                    mod=mod,
                    sim.list=sim.list,
                    group.id=group.id
  )
  
  # Merge model draw outputs into an 3d array
  # This array will have each slice contain a matix with columns for group
  # and rows as posterior samples
  out.array <- abind::abind(out.list,along=3)
  
  # How many groupings exist?
  group.size = dim(out.array)[2]
  
  # Create and merge in group ids to array
  id <- matrix(rep(seq(1,group.size),n.iter),
               nrow = n.iter,
               ncol = group.size,
               byrow = T)
  out.array.ids <- abind::abind(id,out.array,along=3)
  
  # Array needs to be rotated so that each slice is a matrix with group-specific
  # values for each parameter
  out.array.perm <- aperm(out.array.ids,c(2,3,1))
  colnames(out.array.perm) <- c("group_id","Linf","g","inf")
  out.array.perm
}

.single_param_extract <- function(mod,params,sim.list,group.id){
  
  # Check for valid group id
  # if(!group.id %in% c("mu","site","group")) {
  #   stop('group.id must be "mu", "site", or "group"')
  # }
  
  # For parameter grouping or choice, choose the model specific
  # parameter name (e.g., slope for gompertz is g2, and g1 for Von Bert)
  param.mod <-stringr::str_extract(mod@model_pars, "[^_]+$") # (e.g., mu_g1 -> g1)
  param.select <-params[params %in% param.mod]
  
  # Create id for group-specific parameter of interest
  # param.group <-paste(group.id,param.select,sep="_")
  param.group <-paste0(group.id,param.select)
  
  # Extract the posterior draws for each
  out <- sim.list[[param.group]]
  
  # transform to matrix if ony one group (e.g. mu) for consistent formating
  if(length(dim(out)) == 1) out <- matrix(out)
  
  # return parameter matrix
  out
}


# Prediction extracting helpers  -----------------------------------------------

# Extracts posterior distribution of the  predictions of the asymptotic
# and inflection parameters, as well as instantaneous growth predictions given
# specific values of predictors. Predictions are created in Stan model and 
# these functions extracts these predictions into a 3d array (number of 
# predictions, n.parameters X linear predictors, iterations). Each slice is a 
# posterior draw of the prediction.

# REQUIRES: rstan, abind

# extract posterior from model output
.predict_extract <- function(mod.out,mod.dir) {
  
  # Load in model
  mod <- readRDS(file.path(mod.dir,mod.out))
  
  # Select predictions of parameters in actual model
  all_param <- c("pred_Linf","pred_ig")
  param <- all_param[all_param %in% mod@model_pars]
  
  # Extract posteriors
  sim_list <- rstan::extract(mod,param)
  
  # Format each array to bind into one
  pred_list <- purrr::map2(sim_list,names(sim_list), function(x,y) {
    array <- aperm(x,c(2,3,1))  # format array so slices are draws
    colnames(array) <- sapply(seq_len(ncol(array)), function(i) paste0(y,i))  # give each column a unique name
    array
  }
  )
  
  # Bind into one array
  pred_array <-abind::abind(pred_list,along = 2)
  
  # Create a column for prediction id (used for boot summary function)   
  id <- array(rep(seq(1,nrow(pred_array)),dim(pred_array)[3]),
              dim = c(dim(pred_array)[1],1,dim(pred_array)[3]),
              dimnames = list(NULL,"pred_id",NULL)) 
  abind::abind(id,pred_array,along=2)
  
  
}  

# Selects random draws
.linear_predict_sampler <- function (model.out,n.sim) {
  
  # Pull random samples from posterior
  mod_list <- post_draw(model.out,n.sim)
  
  # Extract parameter of interest
  out_list <- lapply(mod_list, function(x) as.data.frame(x))
  abind::abind(out_list,along = 3)
  
}


# Parameter prediction helpers  ------------------------------------------------

# Extracts random posterior draws of group-specific growth parameter
# estimates from n 3d array (groupings) or matrix (no groupings) containing 
# growth parameter estimates from a three parameter growth model. Returns an
# array with group-specific (if present) asymptotic length and
# growth curve inflection parameters, with each slice being a random posterior
# draw.

# REQUIRES: abind

.parameter_sampler <- function (model.out,n.sim,truncate.inf = F,g.mod = NULL) {
  
  # Pull random samples from posterior
  mod_list <- .post_draw(model.out,n.sim)
  
  # Extract parameter of interest
  out_list <- lapply(mod_list, function(x) as.data.frame(x[,c("group_id","Linf","inf")]))
  out <-abind::abind(out_list,along = 3)
  
  # For inf, truncate negative ages to zero, for these fish they 
  # experience greatest growth rate at birth
  if (truncate.inf == F) return(out)
  out[,"inf",][out[,"inf",] <0] <- 0
  out
}


# Curve prediction helpers  ----------------------------------------------------

# Creates group-specific growth curve predictions using length or 
# age values using multiple draws of growth parameter distributions from Stan
# models. Data.frame containing age or length inputs for specified groups are 
# used to create a 3D array containing input data, groupings, and growth 
# predictions. Each slice represents a random posterior draw.

# REQUIRES: dplyr (full), stringr, purrr, abind, parallel (suggested)

.prediction_sampler <-function(model.out,n.sim,...,parallel=F,mc.cores=NULL){
  
  # Parallel processing?
  if (parallel == T) {
    lapply_fun <- function(...) parallel::mclapply(...,mc.cores=mc.cores)
  } else {
    lapply_fun <- lapply
  }
  
  # extract random posterior draws
  mod_list <- .post_draw(model.out,n.sim)
  
  # run specified function using each random growth parameter draws
  out_list <- lapply_fun(mod_list,.growth_predictR,...)
  
  # Bind results into a 3d array
  abind::abind(out_list,along = 3)
  
}

.growth_predictR <- function(model.out,g.mod,input.df,input.var,output.var,wt.df = NULL,dry.wt=1){
  
  # Does model output have groupings?
  if(nrow(model.out)==1){
    # Universal parameter estimates
    out_df <- input.df %>% 
      mutate(Linf=model.out$Linf,
             g = model.out$g,
             inf = model.out$inf)
  } else{
    # Group specific parameters estimates
    out_df <- input.df %>% 
      left_join(model.out,by = join_by(group_id))
  }  
  
  # apply prediction function(s) to data
  preds <- purrr::map(
    output.var,
    ~{
      # Load in helper functions
      # fun_name <- paste(input.var,.x,sep="2")
      fun_name <- paste0(".",input.var,"2",.x)
      fun <- get(fun_name)
      
      # set function arguments
      args <- list(
        input = out_df$input,
        Linf = out_df$Linf,
        g = out_df$g,
        inf = out_df$inf,
        g.mod = g.mod
      )
      
      # add additional arguments for helper functions
      if("days" %in% names(formals(fun))){
        stopifnot("Days missing in input data" = "days" %in% colnames(out_df))
        args$days = out_df$days
      }
      
      if("wt.df" %in% names(formals(fun))){
        stopifnot("Weight-length parameters missing. Please provide dataframe with conversions" = 
                    "days" %in% colnames(out_df))
        args$wt.df <-wt.df
        args$dry.wt <- dry.wt
      }
      
      # apply function
      do.call(fun,args)
    }
  )
  
  # Add prediction columns to data
  names(preds) <- paste0(output.var,"_pred")
  bind_cols(out_df,preds) %>% 
    rename_with(~c(input.var),c(input)) %>% 
    select(-Linf,-g,-inf)
}

# growth prediction helpers  ---------------------------------------------------

# For a given age or length, estimate length, age, instantaneous growth using 
# von Bertalanffy, Gompertz, or logistic growth model parameters.

# REQUIRES: NA

.length2growth <- function(input,g.mod,Linf,g,inf){
  
  stopifnot('Growth model must be from the following models "vb" 
            (von Bertalanffy), "gz" (Gompertz), 
            or "lg" (logistic).'=g.mod %in% c("vb","gz","lg"))  
  
  ## von Bertalanffy  ##
  if(g.mod == "vb"){
    growth = g*(Linf-input)
  }
  
  ## Gompertz  ##
  if(g.mod == "gz"){
    growth <- g*input*log(Linf/input)
  }  # end GZ if statement
  
  ## Logistic  ##
  if(g.mod == "lg"){
    growth = g*input*(1-input/Linf)
  }  # end LG if statement
  
  # Fish above Linf will have negative growth, change this to zero
  growth[growth<0] <- 0
  
  growth
}

.length2age <-function(input,g.mod,Linf,g,inf){
  
  stopifnot('Growth model must be from the following models "vb" 
            (von Bertalanffy), "gz" (Gompertz), 
            or "lg" (logistic).'=g.mod %in% c("vb","gz","lg"))  
  
  ## von Bertalanffy  ##
  if(g.mod == "vb"){
    age = ifelse(input > Linf,
                 Inf,
                 inf-log(1-(input/Linf))/g)
  }
  
  ## Gompertz  ##
  if(g.mod == "gz"){
    age = ifelse(input > Linf,
                 Inf,
                 inf + -log(-log(input/Linf))/g)
  }  # end GZ if statement
  
  if(g.mod == "lg"){
    age = ifelse(input > Linf,
                 Inf,
                 inf-log((Linf/input)-1)/g)
  }  # end LG if statement
  age
}

.age2length <- function(input,g.mod,Linf,g,inf){
  
  stopifnot('Growth model must be from the following models "vb" 
            (von Bertalanffy), "gz" (Gompertz), 
            or "lg" (logistic).'=g.mod %in% c("vb","gz","lg"))  
  
  ## von Bertalanffy  ##
  if(g.mod == "vb"){
    l = Linf * (1 - exp(-g *(input - inf)))
  }
  
  ## Gompertz  ##
  if(g.mod == "gz"){
    l = Linf * exp(-exp(-g * (input - inf)))
  }  # end GZ if statement
  
  ## Logistic  ##
  if(g.mod == "lg"){
    l = Linf/(1 + exp(-g * (input - inf)))
  }  # end LG if statement
  l
}

.age2growth <- function(input,...) {
  l <- .age2length(input,...)
  .length2growth(l,...)
}

# Interval growth helpers  -----------------------------------------------------

# Estimates the exponential growth of a fish during a set interval of time
# given its current length and growth curve outputs. Length is foretasted out
# to the end of the growth interval and Growth rates are provided in 
# terms of fish weight based on provided length-weight parameters

# REQUIRES: NA

.length2interval_growth <- function(input,...,days,wt.df,dry.wt=1){
  
  # Forecast length at end of interval
  length_t <- .length_forcast(input,days,...)
  
  # Estimate growth using weights
  dry_wt = .length2wt(input,wt.df,dry.wt)
  dry_wt_t = .length2wt(length_t,wt.df,dry.wt)
  .exp_growth(dry_wt,dry_wt_t,days)
  
}

.length_forcast <- function(input,days,...) {
  
  # estimate current age
  age = .length2age(input,...)
  
  # Estimate length at end of interval
  # All length >= Linf have inf age, meaning all lengths >=Linf
  # would have length_t = Linf and negative growth. To correct
  # this, change all infinite ages to have length_t=length
  age_t <- age + days
  ifelse(is.infinite(age),
         input,
         .age2length(age_t,...))
}

.length2wt <- function(input, wt.df, dry=1) {
  wt = 10^(wt.df$a + wt.df$b * log10(input*wt.df$c))  # length to weight equation
  dry_wt = wt*dry  # and convert to dry weight (default = 1 so no conversion)
  return(dry_wt)
}

.exp_growth <- function(intial,final,t) {
  log(final/intial)/t
}

# Misc. helpers  ---------------------------------------------------------------

### random posterior draws  ####
# Draw random samples from posterior distribution of growth curve parameters
# from a 3d array (groupings present) or matrix (no groupings). Returns a list
# of growth parameter data.frames (col = parameter, row = groups). Each
# data.frame is a random draw from posterior distribution

.post_draw <- function (model.out,n.sim){
  
  # For outputs without group specific parameters (matrix),
  # each row is a posterior draw
  if(is.na(dim(model.out)[3])){
    
    # How many iterations?
    n.iter = nrow(model.out)
    
    # Random draws
    s <- sample(1:n.iter,n.sim,replace = T)
    
    # Create a list of random draws of growth parameters
    mod_list <- lapply(s,function(x) as.data.frame(model.out[x,]))
  } else {
    
    # For outputs with group specific parameters (arrays),
    # each slice is a posterior draw
    
    # How many iterations?
    n.iter = dim(model.out)[3]
    
    # Random draws
    s <- sample(1:n.iter,n.sim,replace = T)
    
    # Create a list of random draws of growth parameters
    mod_list <- lapply(s, function(x) {
      as.data.frame(matrix(
        model.out[ , , x],
        nrow = dim(model.out)[1],
        dimnames = dimnames(model.out)[1:2]
      ))
    })
  }
  mod_list
}

### Summarize arrays  ###
# Summarizes data in arrays created in bootstrapping or posterior distribution 
# sampling. Returns data.frame containing group specific mean or median values 
# with 95% confidence or credible intervals

# REQUIRES: dplyr (full)

.boot_summary <- function(array,group.var,sum.fun="mean",ci=c(0.0275,0.975)){
  
  # Summary function
  fun <- get(sum.fun)
  
  # summary function label
  fun_label <- paste0("_",sum.fun)
  
  # Calculate summary statistic and 95 CI interval of estimates
  summary <- as.data.frame(apply(array,c(1,2), fun))
  uprs <- as.data.frame(apply(array,c(1,2), quantile, probs=ci[2], na.rm=T))
  lwrs <- as.data.frame(apply(array,c(1,2), quantile, probs=ci[1], na.rm=T))
  
  # Merge into one summary table
  summary %>%
    left_join(lwrs,
              by=group.var,
              suffix = c("","_lwr"))  %>%
    left_join(uprs,
              by=group.var,
              suffix = c(fun_label,"_upr"))
}


