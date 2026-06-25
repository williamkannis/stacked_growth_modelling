#-------------------------------------------------------------------------------
#
#  Growth curve model stacking summary tables
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: 04-13-2026

# DESCRIPTION: Summarizes candidate growth models and stacked growth model 
# outputs, model R2s, loo_cv results, and model stacking weight for results
# tables for manuscript and supporting information.


# House Keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(rstan)
library(parallel)
library(tidyverse)

# Directories
fig_dir <- "figures_cat"
supp_dir <- "figures/supp_tables"
loo_dir <- "loo_outputs_cat"
out_dir <- "stan_outputs/model_out"
input_dir <-"input_data"
fun_dir <- "functions"

# Load in custom functions
source(file.path(fun_dir,"growth_summary_functions.R"))

# Load data
sp_stack_wt <- readRDS(file.path(loo_dir,"stack_wt_out_2026-06-22.rds"))
sp_loo_compare <- readRDS(file.path(loo_dir,"loo_out_2026-06-22.rds"))
r2_df <- readRDS(file.path(loo_dir,"model_r2_2026-06-22.rds"))
stack_param_df <- readRDS(file.path(loo_dir,"stacked_cat_parameters_2026-06-22.rds"))
ind_gmean_df <- readRDS(file.path(loo_dir,"ind_mean_growth_predictions_2026-06-23.rds"))
stack_gmean_df <- readRDS(file.path(loo_dir,"stacked_mean_growth_predictions_2026-06-22.rds"))
age_df <- readRDS(file.path(input_dir,"fsage_cleaned_2026-06-18.rds"))

# Species specific directories
sp <- names(sp_stack_wt)
sp_dir <- sapply(sp,function(x) file.path(out_dir,x))


# Parameter and slope estimates ------------------------------------------------

### Candidate model summary ###

# Extract means and CI
ind_mean_list <- lapply(
  1:length(sp_stack_wt), 
  function (i) {
    mean_ci_batch(
      sp_stack_wt[[i]],
      params = c("cat","sigma_length"),
      sp_dir[i],
      digits=3,
      mc.cores = 6)
    } 
  )
names(ind_mean_list) <- names(sp_stack_wt)

# Add species names
ind_mean_list <- purrr::map2(ind_mean_list,names(ind_mean_list), 
                             function(x,y) x %>% mutate(species =y))

# Bind into one data frame
ind_mean_df <- bind_rows(ind_mean_list)

### Model stacking parameters  ###

# format column names
stack_mean_for <- stack_param_df %>% 
  rename(t_mean = inf_median,
         t_lwr = inf_lwr,
         t_upr = inf_upr)

# Format mean and CIs
stack_mean_list <-lapply(c("Linf","t"),function (i) {
  
  # Extract cat specific parameters
  n_cat <- n_distinct(stack_mean_for$cat_id)
  
  cat_list <-lapply(1:n_cat, function(j){
    
    # Subset grouping
    df_sub <- stack_mean_for %>% 
      filter(cat_id == j)
    
    # Extract species
    sp <- df_sub %>% 
      pull(species)
    
    # Select columns
    df <- df_sub %>% 
      select(-cat_id) %>% 
      select(contains(i))
    
    # ROund digits
    df <-round(df,1)
    
    # Format mean and CI
    mean_ci <- apply(df,1,function(x) paste0(x[1]," (",x[2],", ",x[3],")"))
    mean_ci_df <- data.frame(species = sp,mean_ci)
    colnames(mean_ci_df) <- c("species", paste("cat",i,j,sep="_"))
    mean_ci_df
      
    
  })
  # Merge groupings
  cat_list %>% 
    purrr::reduce(left_join,by="species")
})

# Merge into data frame
stack_mean_df <-stack_mean_list %>% 
  purrr::reduce(left_join,by="species") %>% 
  mutate(model = "stacked")

### Combine stacked and ind parameter means/cis  ###
combined_mean_df <- bind_rows(stack_mean_df,ind_mean_df)


# Inst growth at mean length  --------------------------------------------------

# Format and combine stacking and ind weights
mean_growth_df <-stack_gmean_df %>% 
  mutate(mod = "stacked") %>% 
  bind_rows(ind_gmean_df) %>% 
  
  # Create mean/95 entry for table
  mutate(
    growth_pred_median = round(growth_pred_median,2),
    growth_pred_lwr = round(growth_pred_lwr,2),
    growth_pred_upr = round(growth_pred_upr,2),
    mean_growth = paste0(growth_pred_median," (",growth_pred_lwr,", ",growth_pred_upr,")")) %>% 
  rename(model = mod) %>% 
  select(species,cat_id,model,mean_growth) %>% 
  pivot_wider(names_from = cat_id,values_from = mean_growth) %>% 
  rename(
    mean_growth_1 = `1`,
    mean_growth_2 = `2`,
    mean_growth_3 = `3`
    )
  


# Format Loo tables  -----------------------------------------------------------

# Add species names
sp_loo_compare <- purrr::map2(sp_loo_compare,names(sp_loo_compare), 
                             function(x,y) as.data.frame(x) %>% mutate(species =y))

# Format loo data and estimate CIs for elpd_diff
loo_compare_df <- bind_rows(sp_loo_compare) %>% 
  mutate(lwr = round(elpd_diff - 1.96*se_diff,1),
         upr = round(elpd_diff + 1.96*se_diff,1),
         elpd_diff = round(elpd_diff,1),
         delta_elpd = paste0(elpd_diff, " (",lwr,", ",upr,")")) %>% 
  select(species,elpd_loo,p_loo,se_p_loo,delta_elpd,elpd_diff,se_diff) %>% 
  tibble::rownames_to_column("model")


# Format stacking weights  -----------------------------------------------------

# merge into table
stack_df <- bind_rows(sp_stack_wt) %>% 
    mutate(
      species = substr(model,16,21),
      stack_wt = round(stack_wt,2)
      )
  


# Format and export summarized outputs (Table XX)  -----------------------------

# Join all results
out_table <- combined_mean_df %>% 
  left_join(mean_growth_df) %>% 
  left_join(r2_df) %>% 
  left_join(loo_compare_df) %>% 
  left_join(stack_df) %>% 
  mutate(model = case_when(
    model != "stacked" ~ substr(model,1,2),
    T ~ model
  )) %>% 
  arrange(species,desc(elpd_diff))

# Pretty up NAs
out_table[is.na(out_table)] <- "-"

# Reorder table
out_order <-c("species","model",
              "cat_Linf_1","cat_Linf_2","cat_Linf_3",
              "cat_g_1","cat_g_2","cat_g_3",
              "cat_t_1","cat_t_2","cat_t_3",
              "mean_growth_1","mean_growth_2","mean_growth_3",
              "sigma_length",
              "adj_r2","delta_elpd","stack_wt")
out_table_export <-out_table[,out_order]

# Export
write.csv(out_table_export,file.path(fig_dir,"cat_model_selection_table.csv"))


# Format and export full model outputs (Table sXX)------------------------------

lapply(1:length(sp_stack_wt), function (i) {
  
  # Summarize output
  out<- supp_table_format(
    mod.df = sp_stack_wt[[i]],
    mod.dir = sp_dir[i]
  )
  
  # Export outputs
  fig_name <- paste0("supp_table_",names(sp_stack_wt)[i],".csv")
  write.csv(out,file.path(supp_dir,fig_name),row.names = F)
})


