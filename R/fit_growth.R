#' Stan growth model batch run
#'
#' @description Runs multiple user-specified Bayesian hierarchical growth models 
#' of differing model forms and effect structures.
#' 
#' @param mod.form Vector containing selected growth model forms ("vb" - von 
#' Bertalanffy, "gz" - Gompertz, "lg" - Logistic). Default is c("vb","gz","lg).
#' @param nu Numeric indicating the degrees of freedom for student's t error
#' distributions of lengths. If zero is selected (default), then length error is 
#' modeled using a normal distributions
#' @param fixed.effect Character ("categorical", "continuous","random") indicating 
#' the type of second level effects present in model. Default is no second-level
#' effects (random).
#' @param sample.groups Vector containing the column names used to designate an
#' individual sampling event (e.g site and date) used for random effects. Used
#' to create a single sampling identifier.
#' @param category Name of column containing categorical predictor variable. 
#' Should contain category as a factor. Only necessary if fixed.effects == 
#' "category" Default is NULL.
#' @param predictors Vector with column names of chosen predictor variables. 
#' Only necessary if fixed.effects == "continuous". Default is NULL
#' @param scale T or F: scale and center predictors? Only necessary if 
#' fixed.effects == "continuous". Default is T.
#' @param linear.predictions T or F. Create predictions of growth parameters and
#' and rates across the range of linear predictor variables? Only possible if 
#' fixed.effects == "continuous". Default is F.
#' @param pred.len Number of predictions to make along range of predictor 
#' variables. Only necessary if linear.predictions == T. Default is 100.
#' @param sp Character containing species name for data filtering.
#' @param age.df Data.frame containing age-length data for a set of species, 
#' dates, and sites. Must have columns for species name (species), individual
#' lengths (length), and ages (age). Additional columns for groupings are 
#' required for random effects. If second level fixed effects of selected, data
#' must include columns for categorical groupings, or continuous predictors.
#' @param len.df Data.frame with fish size structure. Used for average length 
#' for inst.growth comparisons across groupings. Must have columns for species 
#' names (species) and lengths (length). If NULL (default), lengths from 
#' age-length data will be used.
#' @param ... Additional arguments to be passed to stan. See documentation for 
#' stan function in rstan.
#' @param parallel T or F. Run multiple model forms in parallel? Warning, if T,
#' it is recommended to set the stan cores augment to zero.
#' @param mc.cores Number of cores used for parallel processing.
#' 
#' @details Reads in pre-created Stan model scripts based on user-specified 
#' growth model forms and fixed effect structures. Models are then ran with the 
#' option of parallel processing. 
#' 
#' Currently, this function can call Bayesian growth models with the three-
#' parameter von Bertalanffy, Gompertz, and logistic growth forms. Length can
#' be modeled using normally distributed errors (nu = 0) or using the student's
#' T distribution with user specified degrees of freedom (nu > 0). MORE DETAILS
#' 
#' Supported effect structures are random effects only, categorical, or continuous
#' second-level effect predictors. In all models, each growth parameter can 
#' vary based on random grouping (e.g. sampling event). For categorical second
#' level effects, random groupings are classified into higher level groups, which
#' can vary in mean growth parameters. For the continuous second-level effects, 
#' differences in random grouping growth parameters are explained by continuous
#' predictors.
#' 
#' The continuous second-level effect model can produce asymptotic length and 
#' length-standardized instantaneous growth predictions of a user-specified 
#' length along the entire range of range of each predictor variable.
#' 
#' Age-at-length and predictor data supplied as a data.frame are automatically  
#' prepared into a named list for Bayesian hierarchical growth modelling in 
#' Stan. Data is filtered by species, assigned a numeric group id based on 
#' user-specified groupings. Mean size of that species across sites is estimated 
#' using an optional size structure data.frame for length-standardized 
#' instantaneous growth estimates. 
#' 
#' If selected, categorical or continuous second-level predictor variables are 
#' prepared. For categorical data, categories are transformed into a numeric 
#' categorical id. For continuous predictors, specified predictor variables are 
#' subset by group id, and optionally scaled and centered. 
#' 
#' @returns Named list containing a list of Stanfit objects for each model 
#' (model_out), a tibble linking sample id to site information (id_bridge), 
#' and mean length used for growth predictions (mean_length). If 
#' fixed.effect == "categorical", output includes an additional list containing 
#' a data.frame (category_labels) that links category names to cat_ids. If 
#' linear.predictions == T, a matrix containing the non-transformed prediction 
#' input values for each predictor (prediction_labels) is included in output 
#' for the use in plot labels.
#' 
#' @export

fit_growth <- function(
    mod.form = c("vb","gz","lg"),
    nu=0,
    fixed.effect = "random",
    sample.groups,
    category = NULL, 
    predictors = NULL,  
    scale = T,
    linear.predictions = F,  
    pred.len = 100,
    sp,   
    age.df, 
    len.df = NULL,
    ...,
    parallel = F,
    mc.cores=NULL){
  
  
  # Prepare input data
  input_list <- .prep_stan_data(
    sp=sp, 
    age.df=age.df, 
    len.df = len.df,
    nu=nu,
    sample.groups = sample.groups,   
    fixed.effect = fixed.effect,
    category = category, 
    predictors = predictors,  
    scale = scale,
    linear.predictions = linear.predictions,  
    pred.len = pred.len
  )
  data <- input_list$stan_data
  
  # Load in appropiate stan scripts
  mods <- sapply(mod.form,.extract_stan_file,fixed.effect)
  
  
  # Run all Stan models in parallel
  if(parallel){  
    out <- parallel::mclapply(
      mods, 
      rstan::stan,
      data=data,
      ...,
      mc.cores = mc.cores
    )
  } else{
    out <- lapply(
      mods, 
      rstan::stan,
      data=data,
      ...,
    )
  }
  names(out) <- gsub("\\..*", "", basename(mods))
  
  # # Sampler diagnostics
  # diag <- lapply(out,.stan_diag)
  
  # return list of model outputs and input data
  # out_diag_list <- list(out,diag)
  # names(out_diag_list) <- c("model_out","sample_diagnostics")
  # c(out_diag_list,input_list) 
  c(list(model_out = out),input_list)

  
}


# Data prep helper function  ---------------------------------------------------
## CHECK FOR LENGTH NAs
.prep_stan_data <- function(
    sp, 
    age.df, 
    len.df = NULL,
    nu = 0,
    sample.groups,
    fixed.effect = "random",
    category = NULL, 
    predictors = NULL,  
    scale = T,
    linear.predictions = F,  
    pred.len = 100){
  
  # Verify age.df has proper coumns and formating
  if(!"species" %in% names(age.df)) {
    stop("Please make sure age.df has a column called 'species'")
  }
  if(!sp %in% age.df$species) {
    stop("Please provide a species name for sp that is included in age.df")
  }
  if(!"length" %in% names(age.df)) {
    stop("Please make sure age.df has a column called 'length'")
  }
  if(!"age" %in% names(age.df)) {
    stop("Please make sure age.df has a column called 'age'")
  }
  if(any(!sample.groups %in% names(age.df))) {
    stop("Please provide column names for sample.groups that appear in age.df")
  }
  
  # If provided, make sure len.df has proper columns and formatting
  if(!is.null(len.df)){
    if(!"species" %in% names(len.df)) {
      stop("Please make sure len.df has a column called 'species'")
    }
    if(!sp %in% len.df$species) {
      stop("Please provide a species name for sp that is included in len.df")
    }
    if(!"length" %in% names(len.df)) {
      stop("Please make sure len.df has a column called 'length'")
    }
  }
  
  # Verify input arguments are valid
  if(!fixed.effect %in% c("random","continuous","categorical")) {
    stop('fixed effect must be "random","continuous",or "categorical"')
  }
  
  # Filter to species of interest and create numeric sample event ids
  sp_df <- age.df %>% 
    dplyr::filter(species == sp) %>% 
    dplyr::group_by(across(tidyr::all_of(sample.groups))) %>% 
    dplyr::mutate(sample_id = dplyr::cur_group_id()) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(sample_id)
  
  # Create table to bridge species specific sample_ids to years and sites
  sample_id_bridge <- sp_df %>% 
    dplyr::distinct(
      dplyr::across(
        tidyr::all_of(c(sample.groups,"species", "sample_id"))
        )
      )
  
  # Average length, to compare growth rates among groupings
  # If no size structure data.set provided, use lengths from age.df
  if(is.null(len.df)) len.df <- age.df
  length_m <- len.df %>% 
    dplyr::filter(species == sp) %>% 
    dplyr::summarise(n = mean(length,na.rm = T)) %>% 
    dplyr::pull()
  
  # Arrange data into input list for Stan analysis
  input_data <- list(
    N = nrow(sp_df),
    N_SITES = dplyr::n_distinct(sp_df$sample_id),
    LENGTH = sp_df$length,
    AGE =sp_df$age,
    ID = sp_df$sample_id,
    LENGTH_M = length_m,
    NU = nu
  )
  
  # Return data if no fixed effects are indicated
  if(fixed.effect == "random") {
    
    if(linear.predictions) {
      stop(
        "Cannot provide linear predictions with the random effect only model"
        )
    }
    
    # Return input data and bridge table as a list
    out <-list(input_data,sample_id_bridge,length_m)
    names(out) <- c("stan_data","id_bridge","mean_length")
    return(out)
  }
  
  # If fixed effects are specified, summarized predictor or categroical data
  # at sampling event level
  sample_df<- sp_df %>% 
    dplyr::select(tidyr::all_of(c("sample_id",category,predictors))) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(sample_id)
  
  # Add categorical second-level predictors if applicable
  if(fixed.effect == "categorical") {
    
    # check if correct type of predictors are provided
    if(!is.null(predictors)) {
    stop("Linear predictors are not possibe with the categorical model")
      }
    if(linear.predictions) {
    stop("Cannot provide linear predictions with the categorical model")
      }
    if(is.null(category)) stop("Please provide category for grouping")

    # select category grouping
    cat_id <- sample_df %>% 
      dplyr::pull(category)
    
    # if category not provided as factor, create factor
    if(!is.factor(cat_id)) cat_id <- factor(cat_id)
    
    # Create numeric id for category groupings
    cat_id <- as.numeric(cat_id)
    
    # link category factor to label
    cat_bridge <- data.frame(
      cat = sample_df[,category],
      cat_id = cat_id
    ) %>% 
      dplyr::distinct()
    
    # Arrange data into input list for Stan analysis
    cat_data <- list(
      N_CAT = dplyr::n_distinct(cat_id),
      CAT = cat_id
      )
    
    # Add categorical data to input list
    input_data <- c(input_data, cat_data)
    
    # Return input data and bridge table as a list
    out <-list(input_data,sample_id_bridge,length_m,cat_bridge)
    names(out) <- c("stan_data","id_bridge","mean_length","category_labels")
    return(out)
  }
  
  # add continuous second-level predictors if applicable
  if(fixed.effect == "continuous"){
    
    # check if correct type of predictors are provided
    if(!is.null(category)) {
    stop("Categorical predictors are not possibe with the continuous model")
      }
    if(is.null(predictors)) stop("Please provide names of continuous predictors")

    # select predictors of choice
    x_df_raw <- sample_df[,predictors]
    
    # Scale and center ?
    if(scale) {
      x_df <- x_df_raw %>% 
        dplyr::mutate(
          dplyr::across(dplyr::everything(),~as.numeric(scale(.x)))
          )
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
      
      # If no predictions are to be made, create empty inputs for stan
      input_data$N_PRED <- 0
      input_data$PRED_X <- matrix(nrow = 0,ncol = input_data$K)
      
      # Return input data and bridge table as a list
      out <-list(input_data,sample_id_bridge,length_m)
      names(out) <- c("stan_data","id_bridge","mean_length")
      return(out)
    }
  }
}


# Stan model look up helper  ---------------------------------------------------

.extract_stan_file <- function(mod, fixed.effects) {
  
  # Create file name
  file <- paste0(mod, "_", fixed.effects, ".stan")
  
  # Load file path from package structure
  file.path <- system.file(
    "stan",
    file,
    # package = utils::packageName()
    package = "growthstack"
  )
  
  if (file.path == "") {
    stop(
      "Stan model not found: ",
      file,
      call. = FALSE
    )
  }
  
  file.path
}

