#-------------------------------------------------------------------------------
# Growth model equation helpers
#-------------------------------------------------------------------------------

# For a given age or length, estimate length, age, instantaneous growth using 
# von Bertalanffy, Gompertz, or logistic growth model parameters.


.length2growth <- function(input,g.mod,Linf,g,inf){
  
  stopifnot('Growth model must be from the following models "vb" 
            (von Bertalanffy), "gz" (Gompertz), 
            or "lg" (logistic).'=g.mod %in% c("vb","gz","lg"))  
  # stopifnot("input length must be postive and non-zero" = all(input > 0))  ## FIX THIS IN TESTING, THIS WARNING IS NOT NECASSARY
  stopifnot("Linf must be postive and non-zero" = Linf > 0)
  stopifnot("g must be postiveand non-zero" = g > 0)
  
  if(any(length(Linf) > 1,length(g) > 1,length(inf) > 1)){
    if(any(
      length(input) != length(Linf),
      length(input) != length(g),
      length(input) != length(inf)
    )){
      stop(
        "if growth parameters provided as vector, must be same length as input"
      )
    }
  }
  
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
  # stopifnot("input length must be postive and non-zero" = all(input > 0))
  stopifnot("Linf must be postive  and non-zero" = Linf > 0)
  stopifnot("g must be postive  and non-zero" = g > 0)
  if(any(length(Linf) > 1,length(g) > 1,length(inf) > 1)){
    if(any(
      length(input) != length(Linf),
      length(input) != length(g),
      length(input) != length(inf)
    )){
      stop(
        "if growth parameters provided as vector, must be same length as input"
      )
    }
  }
  
  # Create index to only estimate valid inputs (those below aysmpote), and 
  # return Inf for those that are invlaid
  age <- rep(Inf,length(input))
  idx <- input <= Linf
  
  if(any(length(Linf) > 1,length(g) > 1,length(inf) > 1)){
    ## von Bertalanffy  ##
    if(g.mod == "vb"){
      age[idx] <- inf[idx]-log(1-(input[idx]/Linf[idx]))/g[idx]
    }
    
    ## Gompertz  ##
    if(g.mod == "gz"){
      age[idx] <- inf[idx] + -log(-log(input[idx]/Linf[idx]))/g[idx]
    }  
    
    if(g.mod == "lg"){
      age[idx] <- inf[idx]-log((Linf[idx]/input[idx])-1)/g[idx]
    }
    age
  }else {
    ## von Bertalanffy  ##
    if(g.mod == "vb"){
      age[idx] <- inf-log(1-(input[idx]/Linf))/g
    }
    
    ## Gompertz  ##
    if(g.mod == "gz"){
      age[idx] <- inf + -log(-log(input[idx]/Linf))/g
    }  
    
    if(g.mod == "lg"){
      age[idx] <- inf-log((Linf/input[idx])-1)/g
    }
  }
  age
}

.age2length <- function(input,g.mod,Linf,g,inf){
  
  stopifnot('Growth model must be from the following models "vb" 
            (von Bertalanffy), "gz" (Gompertz), 
            or "lg" (logistic).'=g.mod %in% c("vb","gz","lg"))  
  # stopifnot("input age must be postive" = all(input >= 0))
  stopifnot("Linf must be postive and non-zero" = Linf > 0)
  stopifnot("g must be postiveand non-zero" = g > 0)
  
  if(any(length(Linf) > 1,length(g) > 1,length(inf) > 1)){
    if(any(
      length(input) != length(Linf),
      length(input) != length(g),
      length(input) != length(inf)
    )){
      stop(
        "if growth parameters provided as vector, must be same length as input"
      )
    }
  }
  
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

.age2interval_growth <- function(input,...,interval,wt.df,dry.wt=1){
  
  # Forecast length at beginning and end of interval
  length <- .age2length(input,...)
  length_t <- .age2length(input+interval,...)
  
  # Estimate growth using weights
  dry_wt = .length2wt(length,wt.df,dry.wt)
  dry_wt_t = .length2wt(length_t,wt.df,dry.wt)
  .exp_growth(dry_wt,dry_wt_t,interval)
}


.length2interval_growth <- function(input,...,interval,wt.df,dry.wt=1){
  
  # Forecast length at end of interval
  length_t <- .length_forcast(input,interval,...)
  
  # Estimate growth using weights
  dry_wt = .length2wt(input,wt.df,dry.wt)
  dry_wt_t = .length2wt(length_t,wt.df,dry.wt)
  .exp_growth(dry_wt,dry_wt_t,interval)
  
}

.length_forcast <- function(input,interval,...) {
  
  ## THIS CAN POSSIBLY BE DISBALED, BUT BLOCK THIS FOR NOW
  stopifnot("Interval must be postive" = interval >=0)
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
  
  stopifnot(
    'wt.df must be data.frame with columns for paremters "a", "b", and "c" ' =
      all(c("a","b","c") %in% names(wt.df))
      )
  stopifnot('wt.df must be a data.frame with 1 row' = nrow(wt.df) == 1)
  stopifnot(
    "dry must be a postive numeric value greater than 0"=
     is.numeric(dry) & dry >0)
  
  wt = 10^(wt.df$a + wt.df$b * log10(input*wt.df$c))  # length to weight equation
  dry_wt = wt*dry  # and convert to dry weight (default = 1 so no conversion)
  return(dry_wt)
}

.exp_growth <- function(intial,final,t) {
  log(final/intial)/t
}

