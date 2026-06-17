#-------------------------------------------------------------------------------
#
#  Stan and loo batch functions
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: Feb 2, 2026

# DESCRIPTION: Contains functions for batch processing of Bayesian hierarchical 
# growth models in Stan. Function prepare data, run models/provide diagnostics,
# and export converged models to a directory. Other functions read in model
# outputs and perform multimodal inference by estimation leave-one-out cross
# validation and estimating stacking weights.

# Requires: dplyr (full), tibble, parallel, rstan (full), loo


# stan_data_prep  --------------------------------------------------------------

#' Stan growth model data preparation
#' 
#' @description Prepares age-length and predictor data for Bayesian hierarchical 
#' growth modelling in Stan.
#' 
#' @param sp Species name for data filtering.
#' @param age.df Data.frame containing age-length data for a set of species, 
#' dates, and sites.
#' @param len.df Data.frame with fish size structure. Used for average length 
#' for inst.growth comparisons across groupings
#' @param fixed.effect Character ("categorical", "linear") indicating the type 
#' of second level effects present in model. Default is no second-level
#' effects (Null).
#' @param pred.df Data.frame with second-level linear or categorical predictors. 
#' Must have matching date and site columns as age.df. Only necessary if 
#' fixed.effects != NULL. Default is NULL
#' @param category Name of column containing categorical predictor variable in
#' pred.df. Should contain catergory as a factor. 
#' Only necessary if fixed.effects == "category" Default is NULL.
#' @param predictors Vector with column names of chosen predictor variables in 
#' pred.df. Only necessary if fixed.effects == "linear". Default is NULL
#' @param scale T or F: scale and center predictors? Only necessary if 
#' fixed.effects == "linear". Default is T.
#' @param linear.predictions T or F. Create predictions of growth parameters and
#' and rates across the range of linear predictor variables? Only possible if 
#' fixed.effects == "linear". Default is F.
#' @param pred.len Number of predictions to make along range of predictor 
#' variables. Only necessary if linear.predictions == T. Default is 100.
#' 
#' @details Data is filtered by species, assigned a numeric group id based 
#' on site and year. Mean size of that species across sites is estimated using
#' an optional size structure data.frame. If selected categorical or linear
#' second-level predictor variables are prepared. For categorical data, categories
#' are transformed into a numeric categorical id. For linear predictors, 
#' predictor variables are subset by group id and selected, and optionally 
#' scaled and centered. The Linear model can produce asymptotic length parameter 
#' and inst. growth of mean size predictions along the entire range of each 
#' predictor variable. This function can prepare input data for predictions 
#' using the full numeric range of provided predictor variables. All data is 
#' then formatted into a name list as per Stan requirements. 
#' 
#' @returns Named list containing the Stan model input data list (stan_data), 
#' a tibble linking sample id to site information (id_bridge), and mean length  
#' used for growth predictions (mean_length). If fixed.effect == "categorical",
#' list contains an additional data.frame (category_labels) that links 
#' category names to cat_ids. If linear.predictions == T, a matrix containing 
#' the non-transformed prediction input values for each predictor 
#' (prediction_labels) is included in output. 
#' 
#' @export

# REQUIRES: dplyr (full)

# MAKE SITE INFO MORE GENERAL

stan_data_prep <- function(sp, age.df, len.df, fixed.effect = NULL, pred.df=NULL, 
                           category = NULL, predictors = NULL, scale = T, 
                           linear.predictions = F,  pred.len = 100){
  
  # Filter to species of interest and create numeric sample event ids
  sp_df <- age.df %>% 
    filter(species == sp) %>% 
    group_by(wateryear,region,site) %>% 
    mutate(sample_id = cur_group_id()) %>% 
    ungroup() %>% 
    arrange(sample_id)
  
  # Create table to bridge species specific sample_ids to years and sites
  sample_id_bridge <- sp_df %>% 
    distinct(wateryear,region,site,species,sample_id)
  
  # Average length, to compare growth rates among groupings
  length_m <- len.df %>% 
    filter(species == sp) %>% 
    summarise(n = mean(length,na.rm = T)) %>% 
    pull()
  
  # Arrange data into input list for Stan analysis
  input_data <- list(
    N = nrow(sp_df),
    N_SITES = n_distinct(sp_df$sample_id),
    LENGTH = sp_df$length,
    AGE =sp_df$ring_count,
    ID = sp_df$sample_id,
    LENGTH_M = length_m
  )
  
  # Return data if no fixed effects are indicated
  if(is.null(fixed.effect)) {
    
    # Return input data and bridge table as a list
    out <-list(input_data,sample_id_bridge,length_m)
    names(out) <- c("stan_data","id_bridge","mean_length")
    return(out)
  }
  
  # Add catergorical second-level predictors if applicable
  if(fixed.effect == "categorical") {
    
    # check if correct type of predictors are provided
    if(is.null(pred.df)) stop("Please provide data.frame containing predictor
                              variables")
    if(!is.null(predictors)) stop("Linear predictors are not possibe with the
                                  categorical model")
    if(linear.predictions) stop("Cannot provide linear predictions with the
                                categorical model")
    if(is.null(category)) stop("Please provide category for grouping")
    
    # Subset predictor data
    sample_df <- sample_id_bridge %>% 
      left_join(pred.df) %>% 
      arrange(sample_id)
    
    # select grouping of choice category and change to factor
    cat_id <- sample_df %>% pull(category)
    
    # if category not provided as factor, create factor
    if(!is.factor(cat_id)) cat_id <- factor(cat_id)
    
    # Create numeric id for category groupings
    cat_id <- as.numeric(cat_id)
    
    # link category factor to label
    cat_bridge <- data.frame(
      cat = sample_df[,category],
      cat_id = cat_id
    ) %>% 
      distinct()
    
    # Arrange data into input list for Stan analysis
    cat_data <- list(
      N_CAT = n_distinct(cat_id),
      CAT = cat_id
      )
    
    # Add categorical data to input list
    input_data <- c(input_data, cat_data)
    
    # Return input data and bridge table as a list
    out <-list(input_data,sample_id_bridge,length_m,cat_bridge)
    names(out) <- c("stan_data","id_bridge","mean_length","category_labels")
    return(out)
  }
  
  # add linear second-level predictors if applicable
  if(fixed.effect == "linear"){
    
    # check if correct type of predictors are provided
    if(is.null(pred.df)) stop("Please provide data.frame containing predictor
                              variables")
    if(!is.null(category)) stop("Categorical predictors are not possibe with 
                                  the linear model")
    if(is.null(predictors)) stop("Please provide names of linear predictors")
    
    
    # Subset predictor data
    sample_df <- sample_id_bridge %>% 
      left_join(pred.df) %>% 
      arrange(sample_id)
    
    # select predictors of choice
    x_df_raw <- sample_df[,predictors]
    
    # Scale and center ?
    if(scale) {
      x_df <- x_df_raw %>% 
        mutate(across(everything(),~as.numeric(scale(.x))))
    } else {
      x_df <- x_df_raw
    }
    
    
    # Arrange data into input list for Stan analysis
    linear_data <- list(
      K = ncol(x_df),
      X = x_df
    )
    
    # Add linear data to input list
    input_data <- c(input_data,linear_data)
    
    # Estimate predictions across linear predictors?
    if (linear.predictions) {
      
      # Prediction input - raw values. Use as labels for pots
      pred_x_df_raw <- apply(
        x_df_raw,
        2,
        simplify = T, 
        function(x) seq = seq(min(x),max(x),length.out = pred.len)
      )
      
      # Scaled and centered prediction input
      pred_x_df <- apply(
        x_df,
        2,
        simplify = T, 
        function(x) seq = seq(min(x),max(x),length.out = pred.len)
      )
      
      # Arrange data into input list for Stan analysis
      pred_data = list(
        N_PRED = nrow(pred_x_df),
        PRED_X = pred_x_df
      )
      
      # Add prediction data to input data
      input_data <- c(input_data,pred_data)
      
      # Return input data and bridge table as a list
      out <-list(input_data,sample_id_bridge,length_m,pred_x_df_raw)
      names(out) <- c("stan_data","id_bridge","mean_length","prediction_labels")
      return(out)
      
    } else {
      # Return input data and bridge table as a list
      out <-list(input_data,sample_id_bridge,length_m)
      names(out) <- c("stan_data","id_bridge","mean_length")
      return(out)
    }
  }
}

# stan_diag_batch --------------------------------------------------------------

#' Batch Stan model run and diagnostics
#'
#' @description Runs multiple Bayesian hierarchical growth models of differing
#' model forms, compiling convergence and sampling diagnostics, and exporting 
#' stanfit objects of converged models
#' 
#' @param mod.forms Vector containing selected  growth model forms ("vb" - von 
#' Bertalanffy, "gz" - Gompertz, "lg" - Logistic). Default is c("vb","gz","lg).
#' @param fixed.effects Types of second level fixed effects ("linear" or 
#' "categorical). Default is NULL, indicating no second-level effects.
#' @param data Named list providing data for stan model. Prepared using 
#' stan_data_prep with the appropriate fixed.effects argument.
#' @param export.dir File path where species specific directory will be created.
#' @param sp name of species used in growth curve. Will be included in Stanfit 
#' file names and the name of new created exported directory.
#' @param mc.core Number of cores used for parallel processing.
#' @param ... Additional arguments to be passed to stan. See documentation for 
#' stan function in rstan.
#' 
#' @details Reads in pre-created Stan model scripts based on user-specified 
#' growth model forms and fixed effect structures. Models are then ran with the 
#' option of parallel processing. Models are then exported to a species 
#' specific directory. Model convergence and sampling diagnostics are provided 
#' in a summarized list.
#' 
#' Currently, this function can call Bayesian growth models with the three-
#' parameter von Bertalanffy, Gompertz, and logistic growth forms. Supported 
#' effect structures are random effects only, categorical, or linear
#' second-level effect predictors. In all models, each growth parameter can 
#' vary based on random grouping (e.g. sampling event). For categorical second
#' level effects, random groupings are classified into higher level groups, which
#' can vary in mean growth parameters. For the linear second-level effects, 
#' differences in random grouping growth parameters are explained by continuous
#' predictors.
#' 
#' Convergence diagnostics include the number of parameters with 
#' rhat > 1.1 and neff > 4*n.chains. Sampler diagnostics include the number of 
#' divergent transitions and tree depths exceeding 10. Loo fit is assessed 
#' using Pareto K values. Pareto K values greater than 1/log(n_eff) or 0.7,
#' which ever is smaller, indicate poor estimation of loo.
#' 
#' @returns List contain convergence and sampling diagnostics for each 
#' input model:model name, number of parameter with rhat < 1.1, number of 
#' parameters with low ess, number of divergent transitions, number of samples 
#' exceeding tree depth, Pareto K values exceeding the threshold, and if the 
#' output was exported.
#' 
#' @export

# REQUIRES: dplyr, parallel,

stan_diag_batch <- function(mod.forms = c("vb","gz","lg"),fixed.effects = NULL,
                            data,sp=NULL,export.dir,...,parallel = F,mc.cores=NULL){
  
  # Model location
  model.dir <- "stan_scripts"
  
  # Do selected model forms and strucutures match avialable choices?
  if(!all(mod.forms %in% c("vb","gz","lg"))) {
    stop('mod.forms must be "vb","gz" or "lg"')}
  if(!is.null(fixed.effects))  {
    if(!fixed.effects %in% c("linear","categorical")){
      stop('fixed.effects must be "linear" or "categorical"')
    }
  }
  
  # List file names based on selected model forms and fixed effect strucutre.
  if(is.null(fixed.effects)) {
    mods <- sapply(mod.forms, function(x) paste0(x,"_random.stan"))
  } else{
    mods <- sapply(mod.forms, function(x) paste0(x,"_",fixed.effects,".stan"))
  }
  
  # Run all Stan models in parallel, export models and creating fit diagnostic 
  # tables
  if(parallel){  
    out <- parallel::mclapply(
      mods, 
      .stan_diag,
      data=data,
      model.dir =model.dir,
      export.dir=export.dir,
      sp=sp,
      ...,
      mc.cores = mc.cores
    )
  } else{
    out <- lapply(
      mods, 
      .stan_diag,
      data=data,
      model.dir =model.dir,
      export.dir=export.dir,
      sp=sp,
      ...,
    )
  }
  
  
  # return list of model diagnostics
  out
  
}


# Helper function
.stan_diag <- function(data,model,model.dir,export.dir,sp,...){
  
  ### Input diagnosis ####
  
  # check that directories exist is in the correct format
  stopifnot("Model directory cannot be found" = dir.exists(model.dir),
            "Export directory cannot be found" = dir.exists(export.dir))
  
  # check if directors have a trailing "/'
  stopifnot('Unwanted "/" at end of model directory' = 
              stringr::str_sub(model.dir,-1,-1) != "/",
            'Unwanted "/" at end of export directory' = 
              stringr::str_sub(export.dir,-1,-1) != "/")
  
  ### Run model ###
  
  # model file.path
  model.full <- file.path(model.dir,model)
  
  # run model
  mod.out <- rstan::stan(
    file = model.full,       # Stan program
    data = data,             # named list of data
    ...,                     # additional arguments (see function documentation)
    refresh = 0             # no progress shown
  )

  ### Check for convergence issues ###
  
  # Extract model summary, including rhats
  mod.summary <- as.data.frame(summary(mod.out)[[1]])
  
  # number of parameters that did not converge (anything >0 is unacceptable)
  no.conv <- nrow(mod.summary %>% filter(Rhat >1.1))
  
  # Bulk effective sample size
  low.eff <- nrow(mod.summary %>% filter(n_eff < 400))
  
  ### Check for sampler issues ###
  
  # Extract sampler parameters after warmup
  sampler.params.post.list <- get_sampler_params(mod.out, inc_warmup = FALSE)
  
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
  loo <- loo::loo(mod.out)
  thres <- sapply(
    1-1/log10((loo$diagnostics$n_eff)), 
    function (x) ifelse(x >.7,.7,x)
    )
  high_k <- loo$diagnostics$pareto_k[loo$diagnostics$pareto_k >thres]

  ### Export model output if convergance is reached ###
  if (no.conv+low.eff +n.diverg == 0){
    
    # Add a species name to file name and directory
    sp_name <- paste0("_",sp)
    export_dir <- file.path(export.dir,sp)
    dir.create(export_dir)
    
    # Create export file names
    export.file <- gsub("\\..*", "", model)
    export_file_sp <- paste0(export.file,sp_name)
    export.dir.file <- file.path(export_dir,export_file_sp)
    
    # Export
    saveRDS(mod.out,paste0(export.dir.file,"_",Sys.Date(),".rds"))
  }
  
  # create a list for convergence and sampling diagnostics
  list(sp = sp,
       model = model,
       high_rhat = no.conv,  # did the model converge
       low_eff = low.eff,  # number of parameters with low ess
       divergence = n.diverg,  # how many divergent transitions were there
       tree_depth = n.tree.depth.10,  # number of samples exceeding tree depth
       high_k = high_k,  # Pareto K values that exceed the threshold
       exported = ifelse(no.conv+low.eff+n.diverg == 0,T,F)  # was a file exported?
       
  )
  
}


# LOO loo_batch  ---------------------------------------------------------------

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

# REQUIRES: parallel loo


loo_batch <- function(out.files = NULL ,out.dir,mc.cores =1) {
  # If only model directory provided, load all models
  if (is.null(out.files)){
    out.files <- list.files(model.dir)
  }
  
  loo_list <- parallel::mclapply(out.files,.loo_import,out.dir)
  names(loo_list) <- out.files
  loo_list
}

# Helper function
.loo_import <- function(out.file,out.dir,cores=1){
  
  # Read model output
  mod_out_path <- file.path(out.dir,out.file)
  out <- readRDS(mod_out_path)
  
  # Extract log likelihood
  log.lik <- loo::extract_log_lik(out,merge_chains = FALSE)
  
  # Calculate loo
  loo::loo(log.lik, cores = cores)
}


# loo_diag  --------------------------------------------------------------------

#' Loo diagnostics
#'
#' @description Runs diagnosis on a list of loo objects.
#' 
#' @param loo_list List containing psis_loo objects
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
  
  # Extract all pareto k values
  loo_len <- length(loo_list)
  k_list <- lapply(seq_len(loo_len),function(x) {
    
    loo <- loo_list[[x]]
    if(is.numeric(k_limit)) thres <- k_limit
    if (k_limit == "ESS") { 
      thres <- sapply(
        1-1/log10((loo$diagnostics$n_eff)), 
        function (x) ifelse(x >.7,.7,x)
        )
      }
    k <-loo$diagnostics$pareto_k
    k[k > thres]
  })
  names(k_list) <- names(loo_list)
  k_list[sapply(k_list,function(x) length(x) > 0)]
  
}


# stack_format  ----------------------------------------------------------------

#' Model stacking formatting
#'
#' @description Estimates model stacking weights, formats results into a 
#' data.frame, and calculates cumulative model weights.
#' 
#' @param loo_list List of psis_loo objects
#' @param cores Number of cores for parallel processing
#' 
#' @details Wrapper for loo::loo_model weights.
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
    rename(stack_wt = x) %>% 
    arrange(desc(stack_wt)) 
  
  # Calculate cumulative model weights
  cum_wt = 0
  for (i in 1:nrow(stack_df)){
    cum_wt <- cum_wt + stack_df$stack_wt[i]
    stack_df[i,"cum_wt"] <- cum_wt
  }
  stack_df
}

