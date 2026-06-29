#-------------------------------------------------------------------------------
# General helper functions
#-------------------------------------------------------------------------------


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


# Summarizes data in arrays created in bootstrapping or posterior distribution 
# sampling. Returns data.frame containing group specific mean or median values 
# with 95% confidence or credible intervals


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
    dplyr::left_join(
      lwrs,
      by=group.var,
      suffix = c("","_lwr")
    )  %>%
    dplyr::left_join(
      uprs,
      by=group.var,
      suffix = c(fun_label,"_upr")
    )
}


