#-------------------------------------------------------------------------------
#
#  Age-class growth curve mode fit
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: 07-05-2026

# DESCRIPTION: Estimates the rsquared values for each model form using the 
# entire age range, the middle 95%, and the upper/lower 2.5% age ranges. 
# Additionally, raw length residuals are plotted. Rsquare tables and residual 
# plots are used for appendices 5 and 6.


# House Keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(dplyr)
library(purrr)
library(rstan)
library(parallel)
library(tidyverse)
library(ggplot2)
devtools::load_all("~/Documents/work/R packages/growthstack")

# Directories
loo_dir <- "outputs/loo_outputs"
out_dir <- "outputs/stan_outputs"
input_dir <-"input_data"
label_dir <- "figures/_labels"
export_dir <- "figures"

# Load in custom functions
# source(file.path(fun_dir,"growth_prediction_functions.R"))

# Load data
sp_stack_wt <- 
  readRDS(file.path(loo_dir,"_cat-stack_wt_out_2026-06-22.rds"))
age_df <- 
  readRDS(file.path(input_dir,"fsage_cleaned_2026-06-18.rds"))
sample_bridge <- 
  readRDS(file.path(label_dir,"fsgwh_sampleid_bridge_2026-08-21.rds")) 
sp_key <-
  readRDS(file.path(label_dir,"fsgwh_sp_key_2026-08-22.rds"))

sp <- names(sp_stack_wt)
sp_dir <- sapply(sp,function(x) file.path(out_dir,x))


# quantile r2 and residuals  ---------------------------------------------------

# Input data
stack.iter <- 10000
sp_list <- rep(sp,each = 4)
quant_list <- rep(c("lwr","mid","upr","all"),length(sp))

# batch species function
fit_list <- purrr::map2(sp_list,quant_list,function(s,q){
  
  # Subset data
  data <- age_df %>% 
    left_join(
      sample_bridge,
      by = join_by(wateryear, region, site, species)
    ) %>% 
    filter(species == s) %>%  
    select(sample_id,length,age)
  
  if(q == "lwr") data <- data %>% 
      filter(age <= quantile(age,0.025)) 
  if(q == "mid") data <- data %>% 
      filter(age > quantile(age,0.025) & age < quantile(age,0.975)) 
  if(q == "upr") data <- data %>% 
      filter(age >= quantile(age,0.975))
  if(nrow(data) < 1) return(NA)
  
  # Candidate model R2
  ind_r2 <- len_R2(
    stack.df = sp_stack_wt[[s]],
    mod.dir = sp_dir[[s]],
    residuals = T,
    data = data,
    stack = F,
    sim = stack.iter,
    sum.fun = "median",
    parallel =T,
    mc.cores = 10
  )
  
  # Stack model R2
  stack_r2 <- len_R2(
    stack.df = sp_stack_wt[[s]],
    mod.dir = sp_dir[[s]],
    residuals = T,
    data = data,
    stack = T,
    sim = stack.iter,
    sum.fun = "median",
    parallel =T,
    mc.cores = 10
  )
  
  # bind into one data.frame
  purrr::map2(
    ind_r2,
    stack_r2,
    function(x,y) rbind(x,y) %>% 
      mutate(species = s,quantile =q)
  )
  
})

# separate rsquared and residual outputs
fit_list_t <- transpose(fit_list)
r2_list <- fit_list_t$rsquared
resid_list <- fit_list_t$residual


# Rsqaured tables (Table s2.1 and s6.2)  ---------------------------------------

# For rsquare table for various summary tables
r2_df <- bind_rows(r2_list[!is.na(r2_list)]) %>% 
  select(-r2) %>% 
  pivot_wider(
    names_from = quantile,
    values_from = adj_r2
  )

saveRDS(
  r2_df,
  file.path(
    loo_dir,
    paste0("_cat-model_r2_",Sys.Date(),".rds")
  )
)

# Export r2 table for model fit appendix results table
r2_table <- r2_df %>% 
  left_join(sp_key) %>% 
  mutate(
    model = casefold(substr(model,1,2),T),
    species = sci_name) %>% 
  select(
    species,
    model,
    lwr,
    mid,
    upr,
    all
  )

write.csv(
  r2_table,
  file.path(
    export_dir,
    "table_s6.2",
    "table_s6.2.csv"
  ),
  row.names = F
)


# Residual plots  --------------------------------------------------------------

resid_df <- bind_rows(resid_list) %>% 
  filter(quantile == "all") %>% 
  mutate(mod = substr(mod,1,2))

lapply(sp, function(s){
  plot_df <- resid_df %>% 
    filter(species ==s)
  
  q1 <- quantile(plot_df$age,.025)
  q3 <- quantile(plot_df$age,.975)
  
  plot <- ggplot(
    plot_df,
    aes(x=age,y=resid,colour = mod)
    )+
    theme_classic(base_size = 40)+
    geom_point(size = 2.5,stroke = 0, shape = 16)+
    geom_hline(yintercept = 0,colour = "red")+
    geom_vline(xintercept = q1,linetype = "dashed")+
    geom_vline(xintercept =q3,linetype = "dashed")+
    xlab("")+
    ylab("") +
    theme(
      legend.position="none",
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 4)
    )
  print(plot)
  
  plot_name <- paste0("_resid_cat_",s,".png")
  ggsave(
    filename = file.path(
      export_dir,
      "figure_s5.2",
      plot_name
      ),
    plot = plot,
    width = 15,
    height = 5,
    dpi = 300
  )
})    

