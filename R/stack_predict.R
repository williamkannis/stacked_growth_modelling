#' Stacked and individual model predictions and growth parameters
#' 
#' @description Create group-specific growth curve predictions 
#' or growth parameters estimates with credible intervals for each candidate
#' model, or a stacked model based on imputed model weights.
#' 
#' @param stack.df Data.frame containing model file names and stacking weights. 
#' Must have columns "model" and "stack_wt. Models must include predictions of
#' instantaneous growth (inst_growth) and asymptotic length (Linf) across a 
#' range of predictor variables.
#' @param mod.dir File path for Stan model output files
#' @param sim Number of posterior draws for stacking
#' @param summarize T or F. Summarize posterior distributions of predictions? If
#' T, returns data.frame. If false, returns 3d array with slice for each 
#' posterior draw. Default is T.
#' @param sum.fun Character ("mean" or "median) for type of summary 
#' statistic of posterior distribution. Default is mean. 
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


stack_predict <- function(
    stack.df, 
    mod.dir,
    sim,
    summarize=T,
    sum.fun=NULL,
    type, 
    group.id,
    pred.input=NULL,
    create.input = T, 
    pred.group = NULL,
    pred.interval =NULL,
    stack=F,
    ...
    ){
  
  ### Prepare inputs ###
  
  # Prediction or parameters?
  fun_name = paste0(".",type,"_sampler")
  .fun = get(fun_name)
  
  # If model stacking, only retain models with stacking weight and sample
  # distributions based on weight.
  if(stack) {
    stack_df <- stack.df %>% 
      dplyr::mutate(n_sim = round(sim*stack_wt)) %>% 
      dplyr::filter(n_sim >0)
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
    ) %>% 
      dplyr::mutate(mod = "stacked")
    
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
      function(x,y) x %>% dplyr::mutate(mod =y)
    )
    
    # Combine into one data.frame
    out_summary <- dplyr::bind_rows(out)
    
  }
  
  # remove temporary prediction id, if applicable
  if(!is.null(out_summary$input_id)){
    out_summary <- out_summary %>% 
      dplyr::select(-input_id)
  } 
  
  # Rename group_id to appropriate grouping name
  group_name <- group.id
  if(group.id == "site") group_name <- "sample_id"
  if(group.id == "cat") group_name <- "cat_id"
  old_name <- "group_id"
  
  out_summary %>% 
    dplyr::rename({{ group_name }} := !! sym(old_name)) 
}




# Model parameter extraction helpers -------------------------------------------

# Extracts posterior distribution of the asymptotic, scaling, and inflection
# parameters from a growth model stanfit object based on user selected
# groupings (grouping - groupt, sampling event - site, or population - mu).
# Returns a 3d array (num groups, 4, iterations). Each slice is a posterior
# draw of the asymptotic, scaling, and inflection parameter values for each
#  group


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
  
  # Extract the posterior draws for each
  out <- sim.list[[param.group]]
  
  # transform to matrix if ony one group (e.g. mu) for consistent formatting
  if(length(dim(out)) == 1) out <- matrix(out)
  
  # return parameter matrix
  out
}




# Parameter prediction helpers  ------------------------------------------------

# Extracts random posterior draws of group-specific growth parameter
# estimates from n 3d array (groupings) or matrix (no groupings) containing 
# growth parameter estimates from a three parameter growth model. Returns an
# array with group-specific (if present) asymptotic length and
# growth curve inflection parameters, with each slice being a random posterior
# draw.


.parameter_sampler <- function (
    model.out,
    n.sim,
    truncate.inf = F,
    g.mod = NULL
    ) {
  
  # Pull random samples from posterior
  mod_list <- .post_draw(model.out,n.sim)
  
  # Extract parameter of interest
  out_list <- lapply(
    mod_list, 
    function(x) as.data.frame(x[,c("group_id","Linf","inf")])
    )
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

.pred_data_prep <- function(
    create.input,
    pred.input,
    max_group,
    pred.group=NULL,
    pred.interval=NULL
    ){
  
  ### Create prediction input data using  user supplied groupings/intervals ###
  if(create.input){
    
    # Create input data for each possible grouping, without suppling interval
    # lengths
    if(all(is.null(pred.group),is.null(pred.interval))) {
      pred_input <- data.frame(group_id = 1:max_group) %>% 
        tidyr::crossing(input = pred.input) %>% 
        dplyr::select(group_id,input)
      pred_input$input_id <- 1:nrow(pred_input)
    }
    
    # Create input data for user supplied groupings without supplying pred.interval
    # lengths
    if(!is.null(pred.group) & is.null(pred.interval)) {
      pred_input <- data.frame(group_id = pred.group) %>% 
        tidyr::crossing(input = pred.input) %>% 
        dplyr::select(group_id,input)
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
        interval = pred.interval
      ) %>% 
        tidyr::crossing(input = pred.input) %>% 
        dplyr::select(group_id,interval,input)
      pred_input$input_id <- 1:nrow(pred_input)
    }
    
    # Create input data for each interval length across each grouping and 
    # input
    if(is.null(pred.group) & !is.null(pred.interval)){
      pred_input <- data.frame(group_id = 1:max_group) %>% 
        tidyr::crossing(
          input = pred.input,
          interval =pred.interval) %>% 
        dplyr::select(group_id,interval,input)
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
        interval = pred.interval,
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


.prediction_sampler <-function(
    model.out,
    n.sim,
    ...,
    parallel=F,
    mc.cores=NULL
    ){
  
  # Parallel processing?
  if (parallel == T) {
    lapply_fun <- function(...) parallel::mclapply(...,mc.cores=mc.cores)
  } else {
    lapply_fun <- lapply
  }
  
  # extract random posterior draws
  mod_list <- .post_draw(model.out,n.sim)
  
  # run specified function using each random growth parameter draws
  out_list <- lapply_fun(mod_list,.growth_predict,...)
  
  # Bind results into a 3d array
  abind::abind(out_list,along = 3)
  
}

.growth_predict <- function(
    model.out,
    g.mod,
    input.df,
    input.var,
    output.var,
    wt.df = NULL,
    dry.wt=1
    ){
  
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
      dplyr::mutate(
        Linf=model.out$Linf,
        g = model.out$g,
        inf = model.out$inf
        )
  } else{
    # Group specific parameters estimates
    out_df <- input.df %>% 
      dplyr::left_join(
        model.out,
        by = dplyr::join_by(group_id)
        )
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
        if(!"interval" %in% colnames(out_df)){
          stop("Please provide pred.interval for interval_growth prediction")
        }
        args$interval = out_df$interval
      }
      
      if("wt.df" %in% names(formals(fun))){
        if(is.null(wt.df)){
          stop(paste0(
            "Weight-length parameters missing.Please provide dataframe with ",
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
  dplyr::bind_cols(out_df,preds) %>% 
    dplyr::rename_with(~c(input.var),c(input)) %>% 
    dplyr::select(-Linf,-g,-inf)
}

