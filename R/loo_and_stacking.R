#-------------------------------------------------------------------------------
#' Loo batch estimation
#' 
#' @description Estimates psis-loo for a list of file names containing
#' Stanfit objects.
#' 
#' @param out.files List of file names for Stanfit objects. If NULL, all 
#' Stanfit objects in directory will be imported
#' @param out.dir Directory containing models of interest.
#' @param mc.cores Number of cores for parallel processing.
#' 
#' @returns A list of psis-objects
#' 
#' @export


loo_batch <- function(out.files = NULL ,out.dir,mc.cores =1) {
  
  # If only model directory provided, load all models
  if (is.null(out.files)){
    out.files <- list.files(model.dir)
  }
  
  # Batch estimation of loo for each model
  loo_list <- parallel::mclapply(out.files,.loo_import,out.dir)
  names(loo_list) <- out.files
  loo_list
}

#  Loo Helper function
.loo_import <- function(out.file,out.dir,cores=1){
  
  # Load in model outputs
  mod_out_path <- file.path(out.dir,out.file)
  out <- readRDS(mod_out_path)
  
  # Calculate loo
  log.lik <- loo::extract_log_lik(out,merge_chains = FALSE)
  loo::loo(log.lik, cores = cores)
}


# ------------------------------------------------------------------------------
#' Model stacking formatting
#'
#' @description Estimates model stacking weights, formats results into a 
#' data.frame, and calculates cumulative model weights.
#' 
#' @inheritParams loo_args
#' 
#' @details Wrapper for loo::loo_model_weights.
#' 
#' @returns Data.frame containing model stacking weights and cumulative 
#' stacking weights for each model 
#' 
#' @export

# OUTPUT: 

# REQUIRES: dplyr (full), loo, tibble


stack_format <- function(loo_list,cores) {
  
  # estimate stacking weights
  stack <- loo::loo_model_weights(loo_list,cores=cores)
  
  # Change stack_wt into dataframe
  stack_df <- as.data.frame(stack) %>% 
    tibble::rownames_to_column("model") %>% 
    dplyr::rename(stack_wt = x) %>% 
    dplyr::arrange(desc(stack_wt)) 
  
  # Calculate cumulative model weights
  cum_wt = 0
  for (i in 1:nrow(stack_df)){
    cum_wt <- cum_wt + stack_df$stack_wt[i]
    stack_df[i,"cum_wt"] <- cum_wt
  }
  stack_df
}

