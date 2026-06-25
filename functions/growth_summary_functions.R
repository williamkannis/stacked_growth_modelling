#-------------------------------------------------------------------------------
#
#  Growth summary statistic functions
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: 6/9/2026

# DESCRIPTION: Functions that summarize the results from Stan model outputs for
# the use in tables and figures.


# mean_ci_batch  ---------------------------------------------------------------

#' Batch median and 95% credible interval summary
#'
#' @description Summarizes the median and 95% credible intervals of posterior
#' distribution of parameters from a set of models for use in tables.
#' 
#' @param mod.df Data.frame containing model file names. Must have column 
#' "model".
#' @param mod.dir File path for Stan model output files
#' @param parallel T or F. Use multiple cores. Only should be used on Linux and 
#' MacOS. Default is F.
#' @param mc.cores Number of core for parallel processing if parallel = T. 
#' Default is NULL.
#' @param params Vector containing parameters to summarize c("mu","tau","beta",
#' "sigma_length").
#' @param digits Number of digits to round values. Default is 3 digits.
#' @param ci Vector containing lower and upper percentiles used for credible 
#' intervals. Default is c(0.025,0.975).
#' 
#' @details Posterior distribution for select parameters are extracted from a 
#' list of models. Posterior distributions are summarized into median and 95% 
#' credible interval. For each parameter median and 95CI values are formatted 
#' into a single character string with a set number of digits.
#' 
#' @returns Data.frame with rows for each model and a single column for each
#' user-specified parameter with median and 95% CI as a character in the 
#' following format: "median (lwr,upr)".
#' 
#' @export 

mean_ci_batch <- function(mod.df,mod.dir,parallel = F,mc.cores=NULL,
                          params = c("mu","tau","beta","sigma_length"),
                          digits=3,ci = c(0.025,0.975)) {
  
  # Load in model names
  mods <- mod.df$model
  
  # Estimate mean and ci
  if (parallel) {
    mean_ci_list <-parallel::mclapply(
      mods,
      .mean_ci_fun,
      mod.dir = mod.dir,
      mc.cores =mc.cores,
      params = params,
      digits=digits,
      ci=ci)
  } else {
    mean_ci_list <-lapply(
      mods,
      .mean_ci_fun,
      mod.dir = mod.dir,
      params = params,
      digits=digits,
      ci=ci)
  }
  
  
  # Combine into single data.frame
  bind_rows(mean_ci_list)
  
}

# Helper function 1
.mean_ci_fun <- function(mod.file,mod.dir,params,digits,ci) {
  
  # Load in model output
  mod <- readRDS(file.path(mod.dir,mod.file))
  
  # Select parameters of interest
  all_params <- mod@model_pars
  params_select <- unlist(sapply(
    params,
    function(i) all_params[grep(i,all_params,T)]
    ))
  params_select <- params_select[!grepl("log",params_select)]
  params_select <- params_select[!grepl("site",params_select)]
  
  # Extract posterior distribution of each parameter
  mod_out <- rstan::extract(mod,params_select)
  
  # Estimate mean and ci of each
  mean_ci_list <- lapply(
    params_select, 
    .mean_ci_helper, 
    mod_out = mod_out,
    mod.file=mod.file, 
    digits = digits,
    ci=ci)
  
  # Combine into one data.frame
  mean_ci_df <-mean_ci_list %>% 
    purrr::reduce(left_join, by = "model")
  
  # Format column names to match across data types
  colnames(mean_ci_df) <- gsub("g1|g2|g3","g",colnames(mean_ci_df))
  colnames(mean_ci_df)<- gsub("t0|ti","t",colnames(mean_ci_df))
  
  if ("tau" %in% colnames(mean_ci_df)){
    mean_ci_df <- mean_ci_df %>% 
      rename(tau_1 = tau)
  }
  
  # Return dataframe
  mean_ci_df
}

# Helper function 2
.mean_ci_helper <- function(param,mod_out,mod.file,digits,ci) {
  
  # Set number of digits
  digit_1 <- c("mu_Linf","mu_ti","mu_t0",
               "cat_Linf","cat_ti","cat_t0")
  digit_2 <- c("sigma_length")
  
  if(param %in% digit_1) digits <- 1
  if(param %in% digit_2) digits <- 2
  
  # Extract parameter of interest
  mod_out <- mod_out[[param]]
  
  # How many estimates per parameters
  out_dim <- length(dim(mod_out))
  
  # Empty list to hold parameter mean and ci
  estimate <-list()
  
  # Extract mean and 95 ci of parameter with single value (e.g., mu)
  if(out_dim == 1) {
    estimate[1] <- median(mod_out)
    estimate[2] <- quantile(mod_out, probs = ci[1])
    estimate[3] <- quantile(mod_out, probs = ci[2])
    
  } 
  
  # Extract mean if more than one value per parameter (e.g., beta[1-n])
  if (out_dim ==2) {
    estimate[[1]] <- apply(mod_out,2,median)
    estimate[[2]] <- apply(mod_out,2,quantile, probs = ci[1])
    estimate[[3]] <- apply(mod_out,2,quantile, probs = ci[2])
    
  }
  
  # Round and retain trailing zeros
  round.id <- paste0("%.",digits,"f") 
  estimate <- lapply(estimate, function(x) sprintf(round.id, x))
  
  # Prepare dataframe with outputs
  if(out_dim == 1) {
    mean_ci <- paste0(estimate[[1]]," (",estimate[[2]],", ",estimate[[3]],")")
    out <- data.frame(mod.file,mean_ci)
    colnames(out) <- c("model",param[1])
    
  } 
  if(out_dim == 2) {
    mean_ci <- sapply(purrr::transpose(estimate), function(x){
      paste0(x[1]," (",x[2],", ",x[3],")")
    })
    out<- data.frame(mod.file,matrix(mean_ci,nrow=1))
    out_names <- sapply(1:length(mean_ci),function(i) paste(param[1],i,sep="_"))
    colnames(out) <- c("model",out_names)
  }
  out
}


# supp_table_format  -----------------------------------------------------------

#' Full Stan model output tables
#'
#' @description Summarizes model outputs for a set of models into one table.
#' 
#' @param mod.df Data.frame containing model file names. Must have column
#' named "model".
#' @param mod.dir File path for Stan model output files
#' 
#' @details Exports a csv file with rows for model parameters for each model
#' type. Includes mean, sd, median, and 95% credible intervals of each 
#' parameter's posterior distribution. Also includes effective sample size and
#' Rhat value.
#' 
#' @export 

supp_table_format <- function(mod.df,mod.dir) {
  
  # Load in model names
  mods <- mod.df$model
  
  # Summarize tables
  format_list <- lapply(mods,.supp_table_format_helper,mod.dir=mod.dir)
  
  # Combine into single data.frame
  bind_rows(format_list) %>% 
    arrange(match(mod,c("vb","gz","lg")))
}

.supp_table_format_helper <- function(mod.dir,mod.file) {
  
  # Load in file and extract model summary
  mod <- readRDS(file.path(mod.dir,mod.file))
  mod_summary <- as.data.frame(summary(mod)[[1]])
  
  # model form
  mod_type <- substr(mod.file,1,2)
  
  # Parameters to summarize
  params <- c("mu_Linf","tau[1]","beta_Linf[1]","beta_Linf[2]","beta_Linf[3]",
              "mu_g1","tau[2]","beta_g1[1]","beta_g1[2]","beta_g1[3]",
              "mu_t0","tau[3]","beta_t0[1]","beta_t0[2]","beta_t0[3]",
              "cor_mat[1,2]","cor_mat[1,3]","cor_mat[2,3]",
              "sigma_length")
  
  # Different model forms have different parameter name variations, 
  # standardize these
  if(mod_type == "gz"){
    params <- gsub("g1","g2",params)
    params <- gsub("t0","ti",params)
  }
  
  if(mod_type == "lg"){
    params <- gsub("g1","g3",params)
    params <- gsub("t0","ti",params)
  }
  
  # Format for supp info table
  mod_summary %>% 
    select(mean,sd,`50%`,`2.5%`,`97.5%`,n_eff,Rhat) %>% 
    tibble::rownames_to_column("parameter") %>% 
    rename(
      median = `50%`,
      lwr = `2.5%`,
      upr = `97.5%`
    ) %>% 
    filter(
      parameter %in% params
    ) %>% 
    arrange(match(parameter,params)) %>% 
    mutate(
      across(-parameter, round, digits = 3),
      mod = mod_type,
    )
}


# beta_mean_ci_batch  ----------------------------------------------------------

#' Batch beta coefficient summary
#'
#' @description Summarizes the median and 95% credible intervals of posterior
#' distribution of beta coefficients from a set of models for use in plotting.
#' 
#' @param stack.df Data.frame containing model file names and stacking weights. 
#' Must have columns "model" and "stack_wt
#' @param wt.cutoff T or F: only include models that have non-zero stacking 
#' weigths? Default is T.
#' @param mod.dir File path for Stan model output files
#' @param ci Vector containing lower and upper percentiles used for credible 
#' intervals. Default is c(0.025,0.975).
#' @param parallel T or F. Use multiple cores. Only should be used on Linux and 
#' MacOS. Default is F.
#' @param mc.cores Number of core for parallel processing if parallel = T. 
#' Default is NULL.
#' 
#' @details Posterior distribution for beta coefficients are extracted from a 
#' list of models. Posterior distributions are summarized into median and 95% 
#' credible interval. 
#' 
#' @returns Data.frame with rows for each model and beta parameter. Has columns
#' mean, lower and upper credible intervals, and model name.
#' 
#' @export 

beta_mean_ci_batch <- function(stack.df, wt.cutoff = T, mod.dir,  
                               ci=c(0.025,0.975),parallel = F,mc.cores=NULL) {
  
  # Filter based on wt
  if(wt.cutoff){
    stack.df <-stack.df %>% 
      mutate(n_samp = stack_wt*1000) %>% 
      filter(n_samp >1)
  }
  
  # Load in model names
  mods <- stack.df$model
  
  # Estimate mean and ci
  if (parallel) {
    mean_ci_list <-parallel::mclapply(
      mods,
      .beta_mean_ci_fun,
      mod.dir = mod.dir,
      mc.cores =mc.cores,
      ci = ci
      )
  } else {
    mean_ci_list <-lapply(
      mods,
      .beta_mean_ci_fun,
      mod.dir = mod.dir,
      ci = ci
        )
  }
  
  # Combine into single data.frame
  bind_rows(mean_ci_list)
  
}

# Helper 1
.beta_mean_ci_fun <- function(mod.file,mod.dir,ci) {
  
  # Load in model output
  mod <- readRDS(file.path(mod.dir,mod.file))
  
  # Select parameters of interest
  all_params <- mod@model_pars
  params_select <- unlist(sapply(
    "beta",
    function(i) all_params[grep(i,all_params,T)])
    )
  params_select <- params_select[!grepl("log",params_select)]
  params_select <- params_select[!grepl("site",params_select)]
  
  # Extract posterior distribution of each parameter
  mod_out <- rstan::extract(mod,params_select)
  
  # Estimate mean and ci of each
  mean_ci_list <- lapply(
    params_select, 
    .beta_helper, 
    mod_out = mod_out,
    mod.file=mod.file, 
    ci = ci
    )
  
  # Combine into one data.frame
  mean_ci_df <-bind_rows(mean_ci_list)
  
  # Format parameter names to match across model types
  mean_ci_df$parameter <- gsub("g1|g2|g3","g",mean_ci_df$parameter)
  mean_ci_df$parameter <- gsub("t0|ti","t",mean_ci_df$parameter)
  
  # Return data.frame
  mean_ci_df
}

# Helper 2
.beta_helper <- function(param,mod_out,mod.file,ci = c(0.025,0.975)) {
  
  # Extract parameter of interest
  mod_out <- mod_out[[param]]
  
  # How many estimates per parameters
  out_dim <- length(dim(mod_out))
  
  # Empty list to hold parameter mean and ci
  estimate <-list()
  
  # Extract mean and 95 ci of parameter with single value (e.g., mu)
  if(out_dim == 1) {
    estimate[1] <- paste0(param[1],"_1")
    estimate[2] <- median(mod_out)
    estimate[3] <- quantile(mod_out, probs = ci[1])
    estimate[4] <- quantile(mod_out, probs = ci[2])
    
  } 
  
  # Extract mean if more than one value per parameter (e.g., beta[1-n])
  if (out_dim ==2) {
    estimate[[1]] <- sapply(1:ncol(mod_out), function(x) paste(param[1],x,sep="_"))
    estimate[[2]] <- apply(mod_out,2,median)
    estimate[[3]] <- apply(mod_out,2,quantile, probs = ci[1])
    estimate[[4]] <- apply(mod_out,2,quantile, probs = ci[2])
    out <- as.data.frame(do.call(cbind,estimate))
    colnames(out) <- c("parameter","mean","lwr","upr")
  }
  out[,"mod_file"] <- mod.file
  out
}
