# len_R2  ----------------------------------------------------------------------

#' R-squared estimation for length-at-age predictions
#' 
#' @description Estimates r-squared of length predictions of all candidate or
#' stacked models. 
#' 
#' @inheritParams stack_predict
#' @param residuals T or F. Return data.frame containing prediction residuals?
#' @param data Data.frame containing actual length-at-age data with appropriate
#' sampling id. Must have columns: age, length, and sample_id.
#' @param stack T or F. Create model stacked predictions or parameter estimates
#' (T), or candidate model specific outputs (F). Default is FALSE
#' @param ... Additional augments passed to [stack_predict].
#'
#' @details Function predicts the length at age for each random grouping. 
#' R-squared is then estimated for each candidate or the stacked model. 
#' 
#' @returns Data.frame containing r-squared and adjusted r-squared for each 
#' model. If 'residuals' == TRUE, returns named list containing r-squared 
#' data.frame, and a data.frame containing predicted values and residuals for 
#' each model.
#' @export

len_R2 <- function(
    stack.df,
    mod.dir,
    sim,
    sum.fun,
    residuals = F,
    data,
    stack,
    ...
    ) {
  
  # model names
  if(stack) mods <- "stacked"
  if(!stack) mods <- stack.df$model
  
  # Predict length at age
  pred_df <-stack_predict(
    stack.df = stack.df, 
    mod.dir = mod.dir,
    type = "prediction",
    group.id = "site",
    pred.input = data$age,
    create.input = F,
    pred.group = data$sample_id,
    stack=stack,
    sim=sim,
    sum.fun = sum.fun,
    input.var = "age",
    output.var = "length",
    ... = ...
    )
  
  # Link predictions to actual data
  linked_df <- pred_df %>% 
    dplyr::distinct(mod,sample_id,age,.keep_all = T) %>% 
    dplyr::left_join(
      data,
      by = dplyr::join_by(sample_id, age),
      relationship = "many-to-many"
    )
  

  # assign mean or median prediction
  if(sum.fun == "mean") linked_df$pred <- linked_df$length_pred_mean
  if(sum.fun == "median") linked_df$pred <- linked_df$length_pred_median
  # if(stack) linked_df$mod <- "stacked"
  
  # R-squared function
  r2_list <- lapply(mods, function (x) {
    
    # subset data for select model
    df <- linked_df %>% dplyr::filter(mod == x)
    # if(stack) df <- linked_df 
    # if(!stack) df <- linked_df %>% dplyr::filter(mod == x)
    
    # # assign mean or median prediction
    # if(sum.fun == "mean") df$pred <- df$length_pred_mean
    # if(sum.fun == "median") df$pred <- df$length_pred_median
    
    # Estimate r2
    reg <- lm(pred ~ length,df)
    data.frame(
      model = x,
      r2 = summary(reg)$r.squared, 
      adj_r2 = summary(reg)$adj.r.squared
    )
  }) 
  r2_df <- dplyr::bind_rows(r2_list)
  
  if(!residuals) return(r2_df)
  
  # Residuals
  res_df <- linked_df %>% 
    dplyr::mutate(resid = length-pred) %>% 
    dplyr::select(mod,age,length,pred,resid)
  
  out <- list(r2_df,res_df)
  names(out) <- c("rsquared","residuals")
  out
  
}
