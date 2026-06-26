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

# curve_predictR  --------------------------------------------------------------

#' Stacked and individual model predictions and growth parameters
#' 
#' @description Create group-specific growth curve predictions 
#' or growth parameters estimates with credible intervals for each candidate
#' model, or a stacked model based on imputed model weights.
#' 
#' 
#' @param stack.df Data.frame containing model file names and stacking wts. 
#' Must have columns "model" and an optional column: "stack_wt if creating 
#' model stacked predictions.
#' @param mod.dir  File path for Stan model output files
#' @param group.id Character ("cat","mu", "site") indicating the grouping 
#' level of model parameters to extract.
#' @param type Character ("parameter" or "prediction"), model stack growth 
#' curve predictions or model parameters parameters
#' @param pred.input vector containing age or length input data if creating
#' predicted growth curves. Default is NULL
#' @param create.input T or F. Create input data in each supplied groupings 
#' and/or interval lengths? If TRUE (default), each user supplied grouping, or 
#' all groupings (if pred.group = NULL) will be assigned each value from 
#' pred.input., If FALSE, pred.input, pred.group, and pred.interval must be 
#' same length.
#' @param pred.group Vector containing numeric identifies for model groupings.
#' Must have equal or less groupings then groupings in model. Groupings 
#' identifiers must match those in the model. If NULL (default), then each model
#' grouping will be assigned each value of pred.input.
#' @param pred.interval Vector containing values for the interval length used
#' to predict interval growth. Only requized if output.vars = "interval_growth".
#' Default is NULL.
#' @param stack T or F. Create model stacked predictions or parameter estimates
#' (T), or candidate model specific outputs (F). Default is FALSE
#' @param sim Number of posterior draws for predictions or parameter estimates
#' @param summarize T or F. Summarize posterior distributions of predictions? If
#' T, returns data.frame. If false, returns 3d array with slice for each 
#' posterior draw. Default is T.
#' @param sum.fun Character ("mean" or "median) for type of summary 
#' statistic of posterior distribution. Default is NULL. 
#' @param ... = Additional arguments passed to .prediction_sampler or 
#' .parameter_sampler auxiliary functions. Required arguments for stacked growth 
#' curve predictions include:
#'   \describe{
#'     \item{input.var}{Character ("length" or "age") to indicate if prediction 
#'     is made using length or age data}
#'     \item{output.var}{Character ("length","age","growth","interval_growth") 
#'     to indicate what type of prediction is being made. If vector is provided, 
#'     multiple prediction columns will be created. Growth indicates 
#'     instantaneous growth, and interval growth is exponential growth during a 
#'     specified interval}
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
#' names provided with optional stacking weights. Parmeters are then either
#' summarized or used to make predicted growth curves using inpute age or 
#' length data. Individual model parameters or predictions then are summarized
#' with the option to combine the posterior distributions into model stacked
#' parameter estimates or growth curves.
#' 
#' @returns Data.frame (Summarized = T) with a row for every input (age or length), 
#' grouping, and model combination (if stack = F). Contains columns for age or length, 
#' model name, grouping index, and prediction summary statistics (mean or 
#' median, lower-lwr, and upper-upr credible interval) for each output variable.
#' Returns a 3 dimensional array (Summarized = F) containing above rows and 
#' columns, with each slice representing one draw from posterior distribution.
#' 
#' @export

# REQUIRES: purrr

curve_predictR <- function(stack.df, mod.dir,type, group.id,pred.input=NULL,
                           create.input = T, pred.group = NULL,
                           pred.interval =NULL,stack=F,sim,
                           summarize=T,sum.fun=NULL, ...){
  
  ### Prepare inputs ###
  
  # Prediction or parameters?
  fun_name = paste0(".",type,"_sampler")
  .fun = get(fun_name)
  
  # If model stacking, only retain models with stacking weight and sample
  # distributions based on weight.
  if(stack) {
    stack_df <- stack.df %>% 
      mutate(n_sim = round(sim*stack_wt)) %>% 
      filter(n_sim >0)
  } else{
    
    # Sampling models evenly
    stack_df <- stack.df %>% 
      mutate(n_sim = sim)
  }
  
  # Model details and number of samples
  mods<- stack_df$model
  sim_list = stack_df$n_sim
  g_mod_list <- substr(mods,1,2)
  
  # Load in parameter samples
  mod_list <- lapply(mods,.param_extract,mod.dir,group.id)
  
  
  ### Compare parameter groupings across models  ###
  
  # number of groupings per model
  n_groups <- sapply(mod_list, function(x) dim(x)[1])
  max_group <- max(n_groups)
  
  # Make sure each model has the same number of groupings if model stacking.
  # Differences in grouping are only acceptable if models that differ only have
  # one grouping (i.e. population mean). These can be duplicated to match the 
  # maximum number of groups
  if(dplyr::n_distinct(n_groups) != 1 & stack){
    
    # How many groups are missing in each model
    group_diff <- max_group -n_groups
    
    # Duplicate population parameters to create equal group sizes for 
    # mode stacking
    mod_list <- purrr::map2(mod_list,group_diff, function(x,y){
      
      # If model has all groups, return model
      if(y==0) return(x)
      
      # How many parameter groups do the models have
      n <- max_group-y
      # If models without grouping have more than one group, hault function
      if(n != 1) stop("Models have unequeal number of groupings")
      
      # Create sequance of missing groups
      group_seq <- (n+1):max_group
      
      # Duplicate parameter value for model by the number of missing groups
      pop_list <-lapply(group_seq, function(z){
        
        # Extract the population parameters
        pop_mat <-x
        
        # change group id to missing groupings
        pop_mat[,"group_id",] <- z
        pop_mat
      })
      abind::abind(c(list(x),pop_list),along=1)
    })
  }
  
  ### Create predictions ###
  
  if(type == "prediction") {
    
    # If groupings are provided, make sure they are groupings in models
    if(!is.null(pred.group)){
      if(max(pred.group) > max_group) {
        stop("User provided groupings exceed number of groupings in model")}
    }
    
    # Create input data based on types of user supplied data
    pred_input <- .pred_data_prep(
      create.input = create.input,
      pred.input = pred.input,
      max_group =max_group,
      pred.group = pred.group,
      pred.interval = pred.interval
    )
    
    # create predictions
    pred_list <-Map(
      function(x,y,z){
        .fun(
          model.out = x,
          g.mod = y, 
          n.sim = z,
          input = pred_input,
          ...=...
        )
      },
      mod_list, 
      g_mod_list, 
      sim_list)
    names(pred_list) <- mods
    
  } else{
    
    ### Extract parameters  ###
    pred_list <-Map(
      function(x,y,z){
        .fun(
          model.out = x,
          ...=...
        )
      },
      mod_list, 
      g_mod_list, 
      sim_list)
    names(pred_list) <- mods
  }
  
  ### Summarize posterior predictive distributions  ###
  
  # Set grouping for summary function
  args <- list(...)
  if("input.var" %in% names(args)){
    group_var <- c("group_id","input_id",args$input.var)
  } else {
    group_var <- c("group_id")
  }
  
  if(stack){
    
    # Bind results into 3d array if stacking
    out <- abind::abind(pred_list, along = 3)
    
    # Return raw values
    if(!summarize) return(out)
    
    # Summarize into data.frame
    out_summary <- .boot_summary(
      out,
      sum.fun=sum.fun,
      group.var=group_var
    )
    
  } else{
    # Return raw values
    if(!summarize) return(pred_list)
    
    # Otherwise summarize models separate and combine into one data.frame
    # Extract means
    out <- lapply(pred_list, .boot_summary,sum.fun=sum.fun, group.var=group_var)
    
    # Add model name
    out <- purrr::map2(
      out,
      names(out),
      function(x,y) x %>% mutate(mod =y)
    )
    
    # Combine into one data.frame
    out_summary <- bind_rows(out)
    
  }
  
  # remove temporary prediction id, if applicable
  if(!is.null(out_summary$input_id)){
    out_summary <- out_summary %>% 
      select(-input_id)
  } 
  
  # Rename group_id to appropriate grouping name
  group_name <- group.id
  if(group.id == "site") group_name <- "sample_id"
  if(group.id == "cat") group_name <- "cat_id"
  old_name <- "group_id"
  
  out_summary %>% 
    rename({{ group_name }} := !! sym(old_name)) 
}


# len_R2  ----------------------------------------------------------------------

#' R-squared estimation for length-at-age predictions
#' 
#' @description Estimates r-squared of length predictions of all candidate or
#' stacked models. 
#' 
#' @param stack.df Data.frame containing model file names and stacking wts. 
#' Must have columns "model" and "stack_wt
#' @param mod.dir File path for Stan model output files
#' @param data Data.frame containing actual length-at-age data with appropriate
#' sampling id. Must have columns: age, length, and sample_id.
#' @param stack T or F. Create model stacked predictions or parameter estimates
#' (T), or candidate model specific outputs (F). Default is FALSE
#' @param sim Number of posterior draws
#' @param sp Character for name of species of which to estimate r-squared
#' @param sum.fun Character ("mean" or "median) for type of summary 
#' statistic of posterior distribution. Default is mean 
#'
#' @details Function predicts the length at age for each random grouping. 
#' R-squared is then estimated for each candidate or the stacked model. 
#' 
#' @returns Data.frame containing r-squared and adjusted r-squared for each 
#' model
#' @export

len_R2 <- function(stack.df,mod.dir,data,stack,sim,sum.fun) {
  
  # model names
  if(stack) mods <- "stacked"
  if(!stack) mods <- stack.df$model
  
  # Predict length at age
  pred_df <-curve_predictR(
    stack.df = stack.df, 
    mod.dir = mod.dir,
    type = "prediction",
    group.id = "site",
    pred.input = data$age,
    pred.group = data$sample_id,
    stack=stack,
    sim=sim,
    sum.fun = sum.fun,
    input.var = "age",
    output.var = "length")

  # Link predictions to actual data
  linked_df <- pred_df %>% 
    left_join(
      data,
      by = join_by(sample_id, age),
      multiple = "first"
    )

  # R-squared function
  r2_list <- lapply(mods, function (x) {
    
    # subset data for select mode
    if(stack) df <- linked_df 
    if(!stack) df <- linked_df %>% filter(mod == x)
    
    # assign mean or median prediction
    if(sum.fun == "mean") df$pred <- df$length_pred_mean
    if(sum.fun == "median") df$pred <- df$length_pred_median
    
    # Estimate r2
    out <- lm(pred ~ length,df)
    data.frame(
      model = x,
      r2 = summary(out)$r.squared, 
      adj_r2 = summary(out)$adj.r.squared
    )
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
  if(!group.id %in% c("mu","site","cat")) {
    stop('group.id must be "mu", "site", or "cat"')
  }

  # For parameter grouping or choice, choose the model specific
  # parameter name (e.g., slope for gompertz is g2, and g1 for Von Bert)
  param.mod <-stringr::str_extract(mod@model_pars, "[^_]+$") # (e.g., mu_g1 -> g1)
  param.select <-params[params %in% param.mod]
  
  # Create id for group-specific parameter of interest
  param.group <-paste(group.id,param.select,sep="_")
  # param.group <-paste0(group.id,param.select)
  
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
  mod_list <- .post_draw(model.out,n.sim)
  
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


# Curve prediction input data helpers  -----------------------------------------

# Takes user inputed arguments and creates a data.frame used with curve 
# prediction helpers to estimate age, length, growth, or interval growth using
# inputed age or length data. Also throws error messages to ensure users have
# supplied the appropriate data for the predictions of their choice.

.pred_data_prep <- function(create.input,pred.input,max_group,pred.group=NULL,pred.interval=NULL){
  
  ### Create prediction input data using  user supplied groupings/intervals ###
  if(create.input){
    
    # Create input data for each possible grouping, without suppling interval
    # lengths
    if(all(is.null(pred.group),is.null(pred.interval))) {
      pred_input <- data.frame(group_id = 1:max_group) %>% 
        tidyr::crossing(input = pred.input) %>% 
        select(group_id,input)
      pred_input$input_id <- 1:nrow(pred_input)
    }
    
    # Create input data for user supplied groupings without supplying pred.interval
    # lengths
    if(!is.null(pred.group) & is.null(pred.interval)) {
      pred_input <- data.frame(group_id = pred.group) %>% 
        tidyr::crossing(input = pred.input) %>% 
        select(group_id,input)
      pred_input$input_id <- 1:nrow(pred_input)
    }
    
    # Create prediction input data using supplied groupings and interval lengths
    if(all(!is.null(pred.group),!is.null(pred.interval))){
      
      # Check if provided data is of the same length
      if(length(pred.group) != length(pred.interval)){
        stop("Please ensure pred.group and pred.interval are the same length")
      }
      
      # Create input data
      pred_input <- data.frame(
        group_id = pred.group,
        pred.interval = pred.interval
      ) %>% 
        tidyr::crossing(input = pred.input) %>% 
        select(group_id,pred.interval,input)
      pred_input$input_id <- 1:nrow(pred_input)
    }
    
    # Create input data for each interval length across each grouping and 
    # input
    if(is.null(pred.group) & !is.null(pred.interval)){
      pred_input <- data.frame(group_id = 1:max_group) %>% 
        tidyr::crossing(
          input = pred.input,
          pred.interval =pred.interval) %>% 
        select(group_id,pred.interval,input)
      pred_input$input_id <- 1:nrow(pred_input)
    }
  } else {
    
    ### Prepare prediction input data using only user supplied data ###
    
    # Check inputs are the same lengths
    if(length(pred.group) != length(pred.input)){
      stop(paste0(
        "pred.group and pred.input must be of the same length if ",
        "create.input = F")
      )
    }
    if(!is.null(pred.interval)){
      if(length(pred.group) != length(pred.interval)){
        stop(paste0(
          "pred.group, pred.input, and pred.interval must all be of the same ",
          "length if create.input = F")
        )
      }
    }
    
    # Prepare data without pred.interval lengths
    if(!is.null(pred.group) & is.null(pred.interval)){
      pred_input <- data.frame(
        group_id = pred.group,
        input = pred.input,
        input_id = 1:length(pred.input)
      )
    }
    
    # Prepare data with pred.interval lengths (for interval growth)
    if(!is.null(pred.group) & !is.null(pred.interval)){
      pred_input <- data.frame(
        group_id = pred.group,
        input = pred.input,
        pred.interval = pred.interval,
        input_id = 1:length(pred.input)
      )
    }
  }
  
  # Return input data.frame
  pred_input
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
  
  # CHeck input data
  if(length(input.var) > 1){
    stop("Please only provide one entry for input.var")
  }
  if(is.null(input.var)){
    stop("Please provide one input variable type")
  }
  if(is.null(output.var)){
    stop("Please provide at least one output variable type")
  }
  
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
      if("interval" %in% names(formals(fun))){
        if("pred.interval" %in% colnames(out_df)){
          stop("Please provide pred.interval  for interval_growth prediction")
        }
        args$interval = out_df$pred.interval
      }
      
      if("wt.df" %in% names(formals(fun))){
        if(is.null(wt.df)){
          stop(paste0(
            "Weight-length parameters missing. Please provide dataframe with ",
            "conversions")
            )
        }
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

.length2interval_growth <- function(input,...,interval,wt.df,dry.wt=1){
  
  # Forecast length at end of interval
  length_t <- .length_forcast(input,interval,...)
  
  # Estimate growth using weights
  dry_wt = .length2wt(input,wt.df,dry.wt)
  dry_wt_t = .length2wt(length_t,wt.df,dry.wt)
  .exp_growth(dry_wt,dry_wt_t,interval)
  
}

.length_forcast <- function(input,interval,...) {
  
  # estimate current age
  age = .length2age(input,...)
  
  # Estimate length at end of interval
  # All length >= Linf have inf age, meaning all lengths >=Linf
  # would have length_t = Linf and negative growth. To correct
  # this, change all infinite ages to have length_t=length
  age_t <- age + interval
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


