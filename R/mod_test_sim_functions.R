#-------------------------------------------------------------------------------
#
#  Simulation and model testing functions
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: 02-27-2026

# DESCRIPTION: Functions used to simulate length-at-age data and to test models
# performance based on simulated data.


# ageLengthSim  ----------------------------------------------------------------

#' Length-at-age data simulation
#'
#' @description Simulates length-at-age and predictor variables using specifed
#' growth forms and effect structures, and user supplied parameter values.
#' 
#' @param sim.input Named list containing the appropriate model parameter values
#' for selected model structure. See details for more information.
#' @param nu Number of degrees of freedom for student's t errors. If zero,
#' error will be simulated using a normal distribution.
#' @param mod.form Vector containing selected growth model forms ("vb" - von 
#' Bertalanffy, "gz" - Gompertz, "lg" - Logistic).
#' @param fixed.effects Character ("categorical", "linear", "none) indicating  
#' the type of second level effects present in model. Default is no second-
#' level effects (none).
#' @param equal.cat T or F. If categorical fixed effects are simulated, should
#' each category be represented by an equal number of sites. Default is T.
#' 
#' @details Users supply a list of growth model parameters to simulate length-
#' at-age data based on a von Bertalanffy, Gompetz, or Logistic growth model.
#' Growth parameters differ based on sampling events using a multivaraite normal
#' distribution, and users are able to specify if growth parameters differ among
#' sampling events based on either linear or categorical predictors. The input 
#' list must contain the following variables:
#' 
#' For all models
#' \describe{
#'     \item{n_sites}{Number of sampling events (i.e., random effect groupings.)}
#'     \item{n_ages}{Number of individuals sampled at each sampling event.}
#'     \item{max_age}{Maximum age of individuals.}
#'     \item{sigma_length}{Length error.}
#'   }
#' 
#' For no fixed effects
#' \describe{
#'     \item{mu_Linf}{Population mean of asymptotic length parameter.}
#'     \item{mu_g}{Population mean of scaling parameter.}
#'     \item{mu_t}{Population mean of inflection parameter.}
#'     \item{tau}{Vector containing the asymptotic length, scaling, and 
#'     infection parameter error.}
#'     \item{cor.Linf.g}{Correlation between asymptotic length and scaling
#'     parameters.}
#'     \item{cor.Linf.t}{Correlation between asymptotic length and inflection
#'     parameters.}
#'     \item{cor.g.t}{Correlation between scaling and inflection parameters.}
#'   }
#'   
#' For linear second-level effects:
#' \describe{
#'     \item{mu_Linf}{Population mean of asymptotic length parameter.}
#'     \item{mu_g}{Population mean of scaling parameter.}
#'     \item{mu_t}{Population mean of inflection parameter.}
#'     \item{tau}{Vector containing the asymptotic length, scaling, and 
#'     infection parameter error.}
#'     \item{cor.Linf.g}{Correlation between asymptotic length and scaling
#'     parameters.}
#'     \item{cor.Linf.t}{Correlation between asymptotic length and inflection
#'     parameters.}
#'     \item{cor.g.t}{Correlation between scaling and inflection parameters.}
#'     \item{beta_Linf}{Vector containing the effect parameters of each 
#'     predictor on asymptotic length. Should be the length of total number of
#'     predictors and match the lengths of beta_t and beta_g.}
#'     \item{beta_t}{Vector containing the effect parameters of each 
#'     predictor on the scaling parameter.}
#'     \item{beta_g}{Vector containing the effect parameters of each 
#'     predictor on the inflection parameter.}
#'   }  
#'   
#' \describe{
#'     \item{cat_Linf}{Population mean of asymptotic length parameter.}
#'     \item{cat_g}{Population mean of scaling parameter.}
#'     \item{cat_t}{Population mean of inflection parameter.}
#'     \item{tau}{Vector containing the asymptotic length, scaling, and 
#'     infection parameter error.}
#'     \item{cor.Linf.g}{Correlation between asymptotic length and scaling
#'     parameters.}
#'     \item{cor.Linf.t}{Correlation between asymptotic length and inflection
#'     parameters.}
#'     \item{cor.g.t}{Correlation between scaling and inflection parameters.}
#'   } 
#'   
#' @returns Named list containing Stan model input data for specified growth
#' form and effect structure. 
#' 
#' @export

#ADD LINEAR PREDICTION INPUTS

ageLengthSim <- function(sim.input,mod.form,nu=0,fixed.effects="none",equal.cat=T) {
  
  ### Prepare input data  ###
  
  # Data structure
  n_sites <- sim.input$n_sites
  n_ages <- sim.input$n_ages
  max_age <- sim.input$max_age
  
  # Log transformed parameters. Leave t0 on normal scale for vb model
  if(fixed.effects != "categorical") {
    mu_log_Linf <- log(sim.input$mu_Linf)
    mu_log_g <- log(sim.input$mu_g)
    if(mod.form != "vb") mu_log_t <- log(sim.input$mu_t)
    if(mod.form == "vb") mu_log_t <- sim.input$mu_t
  }else{
    cat_log_Linf <- log(sim.input$cat_Linf)
    cat_log_g <- log(sim.input$cat_g)
    if(mod.form != "vb") cat_log_t <- log(sim.input$cat_t)
    if(mod.form == "vb") cat_log_t <- sim.input$cat_t
  }
  
  # Mean vector. Estimate means of all categories if using categorical 
  # second level effects 
  if(fixed.effects != "categorical"){
    mu_vector <- c(mu_log_Linf,mu_log_g,mu_log_t)
  }else{
    mu_vector <- c(
      mean(cat_log_Linf),
      mean(cat_log_g),
      mean(cat_log_t)
    )
  }
  
  # variances 
  tau_vector <- sim.input$tau
  
  # Covariance matrix
  cor.mat = matrix(c(1,sim.input$cor.Linf.g,sim.input$cor.Linf.t,
                     sim.input$cor.Linf.g,1,sim.input$cor.g.t,
                     sim.input$cor.Linf.t,sim.input$cor.g.t,1),
                   3, 3)
  cov.mat <-diag(tau_vector) %*% cor.mat  %*%  diag(tau_vector)
  
  ### Simulate random effects ###
  
  # Site random effects
  site_param <-MASS::mvrnorm(n_sites,mu_vector,cov.mat)
  
  # Extract parameter effects
  log_Linf <- site_param[,1]
  log_g <- site_param[,2]
  log_t <- site_param[,3]
  
  ### Simulate second-level fixed effects  ###
  
  # Linear effects
  if(fixed.effects == "linear") {
    
    # number of predictors
    k_sim <- length(sim.input$beta_Linf)
    
    # create random predictor data
    x_sim <- replicate(k_sim, rnorm(n_sites))
    
    # Create site-specific parameter values
    log_Linf <- log_Linf+ x_sim%*%sim.input$beta_Linf
    log_g <- log_g+ x_sim%*%sim.input$beta_g
    log_t <- log_t+ x_sim%*%sim.input$beta_t
  }
  
  # Categorical effects
  if(fixed.effects == "categorical") {
    
    # Number of groupings
    n_cat <- length(cat_log_Linf)
    
    # assign groupings to sites
    # Equal group sizes?
    if(equal.cat){
      cat <- sample(rep(1:n_cat,each = n_sites/n_cat),n_sites)
    } else {
      cat <-sample(1:n_cat,n_sites,replace = T)
    }
    
    # Change groupings into design matrix
    cat_mat <- model.matrix(~ as.factor(cat) - 1)

    # Change group means to group effects
    beta_Linf <- cat_log_Linf - mu_vector[1]
    beta_g <- cat_log_g - mu_vector[2]
    beta_t <- cat_log_t - mu_vector[3]

    # Create site-specific parameter values
    log_Linf <- log_Linf + cat_mat%*%beta_Linf
    log_g <- log_g + cat_mat%*%beta_g
    log_t <- log_t + cat_mat%*%beta_t
  }

  ### Simulate length-at-age data  ###
  
  # Change parameters back to normal scale
  if(mod.form == "vb"){
    param_mat <- cbind(exp(log_Linf),exp(log_g),log_t,1:n_sites)
  }else {
    param_mat <- cbind(exp(log_Linf),exp(log_g),exp(log_t),1:n_sites)
  }
  
  # Create site specific length-age data
  sim_list <-apply(param_mat,1,simplify=F,function(x){
    sample_id = rep(x[4],n_ages)
    age = runif(n_ages,0,max_age)
    
    # Estimate length base on model
    if(mod.form == "vb") length = x[1] * (1 - exp(-x[2] * (age - x[3])))
    if(mod.form == "gz") length = x[1] * exp(-exp(-x[2] * (age - x[3])))
    if(mod.form == "lg") length = x[1]/(1 + exp(-x[2] * (age - x[3])))
    
    # Add error (normal or students t)
    if(nu == 0) length <- length + rnorm(n_ages,0,sim.input$sigma_length)
    if(nu > 0) length <- length + ggdist::rstudent_t(n_ages,nu,0,sim.input$sigma_length)
    cbind(sample_id,length,age)
  })
  
  ### Prepare data for stan input  ###
  
  # Bind simulated lengths into one data.frame
  sim_df <-do.call(rbind,sim_list)
  
  # remove negative lengths
  sim_df <- sim_df[sim_df[,2] > 0,]
  plot(sim_df[,2]~sim_df[,3])
  
  # Create input data
  out <- list(
    N = nrow(sim_df),
    N_SITES = n_distinct(sim_df[,1]),
    LENGTH = sim_df[,2],
    AGE =sim_df[,3],
    ID = sim_df[,1],
    NU=nu,
    LENGTH_M = mean(sim_df[,2])
    )  
  
  # Linear input
  if(fixed.effects == "linear"){
    out$K <- ncol(x_sim)
    out$X <- x_sim
    out$LENGTH_M
    out$N_PRED <- 0
    out$PRED_X <- matrix(nrow=0,ncol=ncol(x_sim))
  }
  
  # Categorical input
  if(fixed.effects == "categorical"){
    out$N_CAT <- n_cat
    out$CAT <- cat
  }
  # return sim data
  out
}

# model_sim_test  --------------------------------------------------------------

#' Growth model simulation testing
#'
#' @description Assess parameter means of credible interval overlap dervied
#' from Stan growth models using simulated data.
#' 
#' @param input Named list providing data for Stan model. Prepared using 
#' ageLengthSim with the appropriate fixed.effects and mod.form argument.
#' @param params Vector containing names of parameter posterior means to 
#' extract. Must match parameter names in Stan model. See details for all
#' possible parameters.
#' @param mod.form Vector containing selected growth model forms ("vb" - von 
#' Bertalanffy, "gz" - Gompertz, "lg" - Logistic).
#' @param fixed.effects Character ("categorical", "linear", "none) indicating  
#' the type of second level effects present in model. Default is no second-
#' level effects (none).
#' @param ... Additional arguments to be passed to stan. See documentation for 
#' stan function in rstan.
#' 
#' @details Calls Bayesian growth models based on user specified model form and
#' second-level effect structures. Models are then fit with simulated data and
#' user-specified parameter means and credible interval overlap with zero are 
#' summarized. If model fails to converge or contains divergent transitions,
#' NA is returned in place of results.
#' 
#' Possible parameters to summarize: "sigma_length", "mu_Linf", "mu_g", "mu_t", 
#' "cat_Linf", "cat_g", "cat_t", "tau", "beta_Linf", "beta_g", and "beta_t".
#' See the details section of ageLengthSim for more parameter descriptions.
#' 
#' @returns Names list containing a table of parameter means ("means") and
#' a table indicating if a  parameter's credible intervals overlap with zero 
#' ("overlap", 1-yes, 0-no).
#' 
#' @export

model_sim_test <- function(input,params,mod.form,fixed.effects="none",...) {
  
  ### Model specific configuration  ###
  if(fixed.effects == "none") fixed.effects <- "random"
  
  # model name
  mod <- paste0(mod.form,"_",fixed.effects,".stan")
  mod_dir <- "stan_scripts"
  
  # growth id
  if(mod.form == "vb") g_id <- 1
  if(mod.form == "gz") g_id <- 2
  if(mod.form == "lg") g_id <- 3
  params[grep("_g", params)] <- paste0(params[grep("_g", params)],g_id)
  
  # inflection id
  if(mod.form == "vb") t_id <- "0" else t_id <- "i"
  params[grep("_t", params)] <- paste0(params[grep("_t", params)],t_id)
  
  ### Run model ###
  out <- stan(
    file = file.path(mod_dir,mod),  # Stan program
    data = input,            # named list of data
    ... = ...
  )
  
  ### Extract model output summary ###
  out_sum <- summary(out)[[1]]
  
  ### Check for sampling issues  ###
  # Check for convergence
  no.converg <- as.data.frame(out_sum) %>% filter(Rhat > 1.1) %>% nrow()
  
  # check for divergent ts
  # Extract sampler parameters after warmup
  sampler.params.post.list <- get_sampler_params(out, inc_warmup = FALSE)
  sampler.params.post.df <-as.data.frame(do.call(rbind,sampler.params.post.list))
  n.diverg <- sum(sampler.params.post.df[,"divergent__"])  # Divergence transitions
  
  # If no convergence do not record values
  if (no.converg+n.diverg > 0) return(NA)
  
  ### Extract values of interest  ###
  mean <-rstan::get_posterior_mean(out,pars = params)[,5]
  
  ### Do CIs overlap with zero?  ###
  overlap_df <- as.data.frame(out_sum) %>% 
    mutate(Overlap0 = case_when(
      `2.5%` < 0 & `97.5%` > 0 ~1,
      T~0
    )) 
  
  overlap <- overlap_df[names(mean),"Overlap0"]
  out <-list(mean,overlap)
  names(out) <- c("mean","overlap")
  out
}

