#-------------------------------------------------------------------------------
#  Stan diagnostic and model selection functions 
#-------------------------------------------------------------------------------

#' Batch Stan model diagnostics
#'
#' @description Runs multiple Bayesian hierarchical growth models of differing
#' model forms, compiling convergence and sampling diagnostics, and exporting 
#' stanfit objects of converged models
#' 
#' @param mod.out stanfit object
#' @inheritParams loo_diag
#' 
#' @details Model convergence and sampling diagnostics are provided 
#' in a summarized list. Convergence diagnostics include the number of 
#' parameters with  rhat > 1.1 and neff > 4*n.chains. Sampler diagnostics 
#' include the number of divergent transitions and tree depths exceeding 10.  
#' Loo fit is assessed using Pareto K values. Pareto K values greater than 
#' 1/log(n_eff) or 0.7, which ever is smaller, indicate poor estimation of loo.
#' 
#' @returns List contain convergence and sampling diagnostics for each 
#' input model:model name, number of parameter with rhat < 1.1, number of 
#' parameters with low ess, number of divergent transitions, number of samples 
#' exceeding tree depth, Pareto K values exceeding the threshold, and if the 
#' output was exported.
#' 
#' @export

stan_diag <- function(mod.out,k_limit){
  
  ### Check for convergence issues ###
  
  # Extract model summary, including rhats
  mod.summary <- as.data.frame(summary(mod.out)[[1]])
  
  # number of parameters that did not converge (anything >0 is unacceptable)
  no.conv <- nrow(mod.summary %>% 
                    dplyr::filter(Rhat >1.1))
  
  # Bulk effective sample size
  low.eff <- nrow(mod.summary %>% 
                    dplyr::filter(n_eff < 400))
  
  ### Check for sampler issues ###
  
  # Extract sampler parameters after warmup
  sampler.params.post.list <- 
    rstan::get_sampler_params(mod.out, inc_warmup = FALSE)
  
  # Convert from list to dataframe
  sampler.params.post.df <-as.data.frame(
    do.call(rbind,sampler.params.post.list)
  )
  
  # Divergence transitions
  n.diverg <- sum(sampler.params.post.df[,"divergent__"])
  
  # Number of tree depths over 10
  n.tree.depth.10  <-nrow(sampler.params.post.df %>%
                            filter(treedepth__>=10))
  
  # Pareto K for loo estimation
  high_k <- .loo_diag_helper(
    loo::loo(mod.out),
    k_limit = k_limit
  )
  
  # create a list for convergence and sampling diagnostics
  list(
    sp = sp,
    high_rhat = no.conv,  # did the model converge
    low_eff = low.eff,  # number of parameters with low ess
    divergence = n.diverg,  # how many divergent transitions were there
    tree_depth = n.tree.depth.10,  # number of samples exceeding tree depth
    high_k = high_k,  # Pareto K values that exceed the threshold
    export = ifelse(no.conv+low.eff+n.diverg == 0,T,F)  # should file be exported?
  )
}


# ------------------------------------------------------------------------------
#' Batch loo diagnostics
#'
#' @description Runs diagnosis on a list of loo objects.
#' 
#' @param loo_list List of psis_loo objects from the loo function.
#' @param k_limit Pareto K thresholds for psis-loo estimation. Default = 0.7.
#' If k_limit = "ESS", threshold will be defined as 1/log(n_eff).
#' 
#' @details Determines if Pareto k statistics are less than a given threshold. 
#' See loo documentation for more information.
#' 
#' @returns Data.frame with model name and max Pareto K for any model beyond 
#' threshold
#' @export

# REQUIRES: NA

loo_diag <- function(loo_list, k_limit = 0.7) {
  
  # Set threshold
  stopifnot('Please provide numeric k_limit or specifiy "ESS".' =
              is.numeric(k_limit) | k_limit == ("ESS"))
  
  # Extract all Pareto k values
  
  k_list <- lapply(loo_list,.loo_diag_helper,k_limit=k_limit)
  loo_len <- length(loo_list)
  # k_list <- lapply(seq_len(loo_len),function(x) {
  #   
  #   loo <- loo_list[[x]]
  #   if(is.numeric(k_limit)) thres <- k_limit
  #   if (k_limit == "ESS") { 
  #     thres <- sapply(
  #       1-1/log10((loo$diagnostics$n_eff)), 
  #       function (x) ifelse(x >.7,.7,x)
  #       )
  #     }
  #   k <-loo$diagnostics$pareto_k
  #   k[k > thres]
  # })
  names(k_list) <- names(loo_list)
  k_list[sapply(k_list,function(x) length(x) > 0)]
  
}


# Loo diagnositc helper
.loo_diag_helper <- function(loo,k_limit) {
  
  # Numeric pareto k thershold?
  if(is.numeric(k_limit)) thres <- k_limit
  
  # THershold based effective sampple size
  if (k_limit == "ESS") { 
    thres <- sapply(
      1-1/log10((loo$diagnostics$n_eff)), 
      function (x) ifelse(x >.7,.7,x)
    )
  }
  
  # Return values above threshold
  k <-loo$diagnostics$pareto_k
  k[k > thres]
}

