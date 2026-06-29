#-------------------------------------------------------------------------------
# Growth model equation helpers
#-------------------------------------------------------------------------------

# For a given age or length, estimate length, age, instantaneous growth using 
# von Bertalanffy, Gompertz, or logistic growth model parameters.


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