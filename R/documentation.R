#' Stan growth model arguments 
#' 
#' @param mod.forms Vector containing selected growth model forms ("vb" - von 
#' Bertalanffy, "gz" - Gompertz, "lg" - Logistic). Default is c("vb","gz","lg).
#' @param NU Numeric indicating the degrees of freedom for student's t error
#' distributions of lengths. If zero is selected (default), then length error is 
#' modeled using a normal distributions
#' @param fixed.effect Character ("categorical", "linear","random") indicating 
#' the type of second level effects present in model. Default is no second-level
#' effects (random).
#' @param sample.groups Vector containing the column names used to designate an
#' individual sampling event (e.g site and date) used for random effects. Used
#' to create a single sampling identifier.
#' @param category Name of column containing categorical predictor variable. 
#' Should contain category as a factor. Only necessary if fixed.effects == 
#' "category" Default is NULL.
#' @param predictors Vector with column names of chosen predictor variables. 
#' Only necessary if fixed.effects == "linear". Default is NULL
#' @param scale T or F: scale and center predictors? Only necessary if 
#' fixed.effects == "linear". Default is T.
#' @param linear.predictions T or F. Create predictions of growth parameters and
#' and rates across the range of linear predictor variables? Only possible if 
#' fixed.effects == "linear". Default is F.
#' @param pred.len Number of predictions to make along range of predictor 
#' variables. Only necessary if linear.predictions == T. Default is 100.
#' @param sp Character containing species name for data filtering.
#' @param age.df Data.frame containing age-length data for a set of species, 
#' dates, and sites. Must have columns for species name (species), individual
#' lengths (length), and ages (age). Additional columns for groupings are 
#' required for random effects. If second level fixed effects of selected, data
#' must include columns for categorical groupings, or linear predictors.
#' @param len.df Data.frame with fish size structure. Used for average length 
#' for inst.growth comparisons across groupings. Must have columns for species 
#' names (species) and lengths (length). If NULL (default), lengths from 
#' age-length data will be used.
#' 
#' @name growth_mod_args
NULL

#' Loo function arguments
#' 
#' @param loo_list List of psis_loo objects from the loo function.
#' @param cores Number of cores for parallel processing in the 
#' loo_model_weights function in loo package.
#' @param k_limit Pareto K thresholds for psis-loo estimation. Default = 0.7.
#' If k_limit = "ESS", threshold will be defined as 1/log(n_eff).
#' 
#' @name loo_args
NULL