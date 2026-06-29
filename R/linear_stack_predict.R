#' Model-stacked instantaneous growth predictions across values of predictors
#'
#' @description Creates model-stacked predictions of instantaneous growth rate
#' across a range of predictor variables.
#' 
#' @inheritParams stack_predict
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


linear_stack_predict <- function(stack.df, mod.dir, sim,sum.fun){
  
  # Retain models with 
  stack <- stack.df %>% 
    dplyr::mutate(n_sim = round(sim*stack_wt)) %>% 
    dplyr::filter(n_sim >0)
  
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
  
# Prediction extracting helpers  -----------------------------------------------

# Extracts posterior distribution of the  predictions of the asymptotic
# and inflection parameters, as well as instantaneous growth predictions given
# specific values of predictors. Predictions are created in Stan model and 
# these functions extracts these predictions into a 3d array (number of 
# predictions, n.parameters X linear predictors, iterations). Each slice is a 
# posterior draw of the prediction.


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
  mod_list <- .post_draw(model.out,n.sim)
  
  # Extract parameter of interest
  out_list <- lapply(mod_list, function(x) as.data.frame(x))
  abind::abind(out_list,along = 3)
  
}