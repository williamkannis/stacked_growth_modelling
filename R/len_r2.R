# len_R2  ----------------------------------------------------------------------

#' R-squared estimation for length-at-age predictions
#' 
#' @description Estimates r-squared of length predictions of all candidate or
#' stacked models. 
#' 
#' @inheritParams curve_predictR
#' @param data Data.frame containing actual length-at-age data with appropriate
#' sampling id. Must have columns: age, length, and sample_id.
#' @param stack T or F. Create model stacked predictions or parameter estimates
#' (T), or candidate model specific outputs (F). Default is FALSE
#'
#' @details Function predicts the length at age for each random grouping. 
#' R-squared is then estimated for each candidate or the stacked model. 
#' 
#' @returns Data.frame containing r-squared and adjusted r-squared for each 
#' model
#' @export

len_R2 <- function(
    stack.df,
    mod.dir,
    sim,
    sum.fun,
    data,
    stack
    ) {
  
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
    dplyr::left_join(
      data,
      by = dplyr::join_by(sample_id, age),
      multiple = "first"
    )
  
  # R-squared function
  r2_list <- lapply(mods, function (x) {
    
    # subset data for select mode
    if(stack) df <- linked_df 
    if(!stack) df <- linked_df %>% dplyr::filter(mod == x)
    
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
  dplyr::bind_rows(r2_list)
}
