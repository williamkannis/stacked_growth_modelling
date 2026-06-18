#-------------------------------------------------------------------------------
#
#  Growth curve model plotting
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: 04-13-2026

# DESCRIPTION: Creates length/growth curve plots, predicted growth at hydrology 
# plots, and model coefficient plots for manuscript.


# House Keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(rstan)
library(parallel)
library(tidyverse)

# Directories
out_dir <- "stan_outputs/model_out"
loo_dir <- "loo_outputs"
pred_dir <- "input_data"
plot_dir <- "stan_outputs/plotting_info"
export_dir <- "figures"
mu_dir <- file.path(export_dir,"mu_plots")
fun_dir <- "functions"

# Custom functions
source(file.path(fun_dir,"growth_summary_functions.R"))

# Load data
sp_stack_wt <- readRDS(file.path(loo_dir,"stack_wt_out_2026-06-16.rds"))
curve_df <- readRDS(file.path(loo_dir,"stacked_curves_2026-06-16.rds"))
mu_curve_df <- readRDS(file.path(loo_dir,"stacked_mu_curves_2026-06-16.rds"))
ind_mu_curve_df <- readRDS(file.path(loo_dir,"ind_mu_curves_2026-06-16.rds"))
pred_bridged <- readRDS(file.path(loo_dir,"stacked_growth_predictions_2026-06-16.rds"))
pred_df <-readRDS(file.path(pred_dir,"fsgrw_predictors_2026-04-23.rds"))
age_df <- readRDS(file.path(pred_dir,"fsage_cleaned_2026-06-18.rds"))

# Species specific directories
sp <- names(sp_stack_wt)
sp_dir <- sapply(sp,function(x) file.path(out_dir,x))


# Prepare data  ----------------------------------------------------------------

# Actual data for points and rug lines
actual_df <- age_df %>% left_join(pred_df) %>% 
  # filter(species != "JORFLO") %>% 
  mutate(group_name = paste(region,site,wateryear)) %>% 
  select(species,group_name,age,length,PC1,PC2,PC3)

# All plot parameters
sp.colors <- c("#b5a331","#339d38","#c26a77","#8c6d3f","#2f2585","#2b695c")
line.size = 1
rug.size = 3
ribbon.alpha = 0.1
base.size = 45
border.size =3
point.size <- 1


# mu Plot parameters
max_age_mu <-plyr::round_any(max(actual_df$age),5,ceiling)
max_age_mu <- 360
mu_breaks <- seq(0,max_age_mu,length.out =4)
base.size.mu <-45


# Population curves - model comparisons (Fig 3)  -------------------------------

# Prepare data
mod_compare_df <- mu_curve_df %>% 
  mutate(mod = "stacked") %>% 
  filter(!species %in% c("POELAT","HETFOR")) %>% 
  bind_rows(ind_mu_curve_df) %>% 
  mutate(mod = substr(mod,1,2),
         mod = factor(mod,levels = c("vb","gz","lg","st")))

# Species specific plotting loop
for(i in 1:length(sp)){
  
  # Subset based on species
  mod_compare_sp <- mod_compare_df %>% 
    filter(species ==sp[i])
  
  # Set colors for top models vs other models
  # default is to highlight stacked models
  mod_compare_sp <- mod_compare_sp %>% 
    mutate(top_mod = case_when(
      mod == "st" ~ 1,
      T~0
    ))
  
  # In the species without stacking, highlight top model
  if (sp[i] == "POELAT") mod_compare_sp$top_mod[mod_compare_sp$mod == "lg"] <- 1
  if (sp[i] == "HETFOR") mod_compare_sp$top_mod[mod_compare_sp$mod == "vb"] <- 1
  
  # Assign color scheme
  mod_compare_sp$top_mod <- factor(mod_compare_sp$top_mod,levels = c(0,1))
  top_mod_color <- c("black","red")
  
  # Assign line types for each model
  # Number of models, if 3 opposed to 4, then no stacked model
  n.mods <- n_distinct(mod_compare_sp$mod)
  line_vec <- c(2,3,4,1)[1:n.mods]
  
  # Length at age plot
  length_compare <-ggplot(data=mod_compare_sp,aes(x = age,y=length_pred_median,group = mod,linetype = mod,colour = top_mod))+
    geom_line() +
    geom_ribbon(aes(ymin=length_pred_lwr, ymax=length_pred_upr,colour = top_mod), linetype=0, alpha=0.1)+
    # geom_point(data = actual_df %>% filter(species == sp[i]),aes(x=age,y = length),size = point.size,inherit.aes = F)+
    geom_rug(data=actual_df %>% filter(species == sp[i]),aes(x=age),inherit.aes  = F)+
    scale_linetype_manual(values =line_vec) +
    scale_color_manual(values =top_mod_color)+
    xlab("")+
    ylab("")+
    xlim(c(0,max_age_mu))+
    theme_classic(base_size = base.size.mu)+ 
    theme(legend.position="none") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = border.size))
  print(length_compare)

  # Growth at age plot
  growth_compare <- ggplot(data=mod_compare_sp,aes(x = age,y=growth_pred_median,group = mod,linetype = mod, colour = top_mod))+
    geom_line() +
    geom_ribbon(aes(ymin=growth_pred_lwr, ymax=growth_pred_upr,group = mod,colour = top_mod), linetype=0, alpha=0.1)+
    geom_rug(data=actual_df %>% filter(species == sp[i]),aes(x=age),inherit.aes  = F)+
    scale_linetype_manual(values =line_vec) +
    scale_color_manual(values =top_mod_color)+
    xlab("")+
    ylab("")+
    xlim(c(0,max_age_mu))+
    theme_classic(base_size = base.size.mu)+ 
    theme(legend.position="none") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = border.size))
  print(growth_compare)
  
  # Export length
  length_compare_name <- paste0("length_age_compare_plot_",sp[i],".png")
  ggplot2::ggsave(file.path(mu_dir,length_compare_name),
                  length_compare,
                  width = 10,
                  height = 8,
                  dpi = 300)

  # Export growth
  growth_compare_name <- paste0("growth_age_compare_plot_",sp[i],".png")
  ggplot2::ggsave(file.path(mu_dir,growth_compare_name),
                  growth_compare,
                  width = 10,
                  height = 8,
                  dpi = 300)
  
}


# Sampling event curves (Fig 4)  -----------------------------------------------

# Set gradient for colors
min_grad <- pred_bridged %>% select(PC1) %>% min()
max_grad <- pred_bridged %>% select(PC1) %>% max()
min_col <- "#132B43"
max_col <- "#56B1F7"

# Set y limits for length plots
age_limit <- actual_df %>% 
  filter(species != "JORFLO") %>% 
  group_by(species) %>% 
  summarize(age=max(age)) 
max_length <- curve_df %>% 
  right_join(age_limit,join_by(species,age)) %>% 
  group_by(species) %>% 
  summarize(length = max(length_pred_upr)) %>% 
  ungroup() %>% 
  summarize(length = max(length)) %>% 
  pull(length) %>% 
  plyr::round_any(5,ceiling)
max_length <- 60

# length breaks
len_break <- seq(0,max_length,length.out = 4)

# Set y limits for growth plots
max_growth <- curve_df %>% 
  summarize(length = max(growth_pred_upr)) %>% 
  ungroup() %>% 
  summarize(length = max(length)) %>% 
  pull(length) %>% 
  plyr::round_any(.1,ceiling)

# Growth breaks
growth_break <- seq(0,max_growth,length.out = 4)
grow_y_lim_adj <-c(0,0,0,.1,.1,0)

# Species specific length plots
for(i in 1:6){
  
  # Subset actual data for ruglines
  l_fish_plot <- actual_df %>% 
    filter(species == sp[i]) 
  
  # Limit plot only predictions within the range of actual fish sampled
  max_age <-plyr::round_any(max(l_fish_plot$age),25,f=ceiling)
  max_age <- 360
  l_plot_df <- curve_df %>% 
    filter(!is.na(sample_id),
           species == sp[i],
           age <=max_age) %>% 
    left_join(l_fish_plot %>% distinct(group_name,PC1))
  
  # Set y axis limits
  max_length <- plyr::round_any(max(l_plot_df$length_pred_upr),5,ceiling)
  
  # Create plot
  p <-ggplot(data=l_plot_df,aes(x = age,y=length_pred_median,group = group_name,colour = PC1))+
    geom_line() +
    geom_ribbon(aes(ymin=length_pred_lwr, ymax=length_pred_upr,fill = PC1), linetype=0, alpha=0.1)+
    geom_point(data = l_fish_plot,aes(x=age,y = length,group = group_name, colour = PC1),size = point.size)+
    scale_color_gradient(
      low = min_col,
      high = max_col,
      limits = c(min_grad, max_grad))+
    scale_fill_gradient(
      low = min_col,
      high = max_col,
      limits = c(min_grad, max_grad))+
    xlab("")+
    ylab("")+
    scale_y_continuous(
      limits = c(0,max_length),
      #breaks = len_break
      )+
    xlim(c(0,max_age))+
    theme_classic(base_size = base.size)+ 
    theme(legend.position="none") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = border.size)) 
  print(p)
  
  # Create file name based on species
  length_plot_name <- paste0("length_at_age_plots/",sp[i],"_length_age_plot.png")
  
  # Export plot
  ggplot2::ggsave(file.path(export_dir,length_plot_name),
                  p,
                  width = 10,
                  height = 8,
                  dpi = 300)
  
}


# Growth at length plot - species specific
for(i in 1:6){
  
  # Actual data
  g_fish_plot <- actual_df %>% 
    filter(species == sp[i]) 
  
  # Limit plot only predictions within the range of actual fish sampled
  max_age <-plyr::round_any(max(g_fish_plot$age),25,f=ceiling)
  max_age <- 360
  g_plot_df <- curve_df %>% 
    filter(!is.na(sample_id),
           species == sp[i],
           age <= max_age)%>% 
    left_join(g_fish_plot %>% distinct(group_name,PC1)) 
  
  # Set y axis limits
  max_growth <- plyr::round_any(max(g_plot_df$growth_pred_upr),.1,ceiling)
  max_growth <- max_growth+grow_y_lim_adj[i]
  #growth_break <- seq(0,max_growth,length.out = 4)
  #growth_break <- round(growth_break,2)

  # Plotting
  library(ggplot2)
  p <-ggplot(data=g_plot_df,aes(x = age,y=growth_pred_median,group = group_name,color=PC1))+
    geom_line() +
    geom_ribbon(aes(ymin=growth_pred_lwr, ymax=growth_pred_upr,group = group_name,fill=PC1), linetype=0, alpha=0.1)+
    geom_rug(data=g_fish_plot,aes(x=age,colour = PC1),inherit.aes  = F,line.width = rug.size)+
    scale_color_gradient(
      low = min_col,
      high = max_col,
      limits = c(min_grad, max_grad))+
    scale_fill_gradient(
      low = min_col,
      high = max_col,
      limits = c(min_grad, max_grad))+
    xlab("")+
    ylab("")+
    xlim(c(0,max_age))+
    scale_y_continuous(
      limits = c(0,max_growth),
      #breaks = growth_break
    )+
    theme_classic(base_size = base.size)+ 
    theme(legend.position="none") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = border.size)) 
  print(p)
  
  # Create file name based on species
  growth_plot_name <- paste0("growth_at_age_plots/",sp[i],"_growth_age_plot.png")
  
  # export
  ggplot2::ggsave(file.path(export_dir,growth_plot_name),
                  p,
                  width = 10,
                  height = 8,
                  dpi = 300)
  
}



# Population curves - all species (Fig 5)  -------------------------------------

# Age at length
mu_l_plot <-ggplot(data=mu_curve_df,aes(x = age,y=length_pred_median,group = species,colour = species))+
  geom_line() +
  geom_ribbon(aes(ymin=length_pred_lwr, ymax=length_pred_upr,fill = species), linetype=0, alpha=0.1)+
  geom_rug(data=actual_df,aes(x=age,colour = species),inherit.aes  = F)+
  # geom_point(data = actual_df,aes(x=age,y = length,group = group_name, colour = species),size = point.size)+
  scale_color_manual(values =sp.colors)+
  xlab("")+
  ylab("")+
  xlim(c(0,max_age_mu))+
  theme_classic(base_size = base.size.mu)+ 
  theme(legend.position="none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = border.size))
ggplot2::ggsave(file.path(mu_dir,"mu_length_age_plot.png"),
                mu_l_plot,
                width = 10,
                height = 8,
                dpi = 300)


# growth at age
mu_g_plot <-ggplot(data=mu_curve_df,aes(x = age,y=growth_pred_median,group = species,color=species))+
  geom_line() +
  geom_ribbon(aes(ymin=growth_pred_lwr, ymax=growth_pred_upr,group = species,fill=species), linetype=0, alpha=0.1)+
  geom_rug(data=actual_df,aes(x=age,colour = species),inherit.aes  = F)+
  scale_color_manual(values =sp.colors)+
  xlab("")+
  ylab("")+
  xlim(c(0,max_age_mu))+
  theme_classic(base_size = base.size.mu)+ 
  theme(legend.position="none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = border.size)) 
ggplot2::ggsave(file.path(mu_dir,"mu_growth_age_plot.png"),
                mu_g_plot,
                width = 10,
                height = 8,
                dpi = 300)


# Beta coefficient plot (Fig 6)  -----------------------------------------------

# Remove JORFLO
sp_stack_wt_for <- sp_stack_wt[names(sp_stack_wt) != "JORFLO"]
sp_dir_for <- sp_dir[names(sp_dir) != "JORFLO"]

# Extract beta parameters and CIs
beta_list <- lapply(1:length(sp_stack_wt_for), function (i) {
  beta_mean_ci_batch(sp_stack_wt_for[[i]],T,sp_dir_for[i],mc.cores = 3)
} )

# Format beta data
beta_df <- bind_rows(beta_list) %>% 
  mutate(
    mean = as.numeric(mean),
    lwr = as.numeric(lwr),
    upr = as.numeric(upr),
    species = substr(mod_file,8,13),
    mod = substr(mod_file,1,2),
    mod = factor(mod,levels = c("vb","gz","lg")),
    response = str_extract(parameter, "(?<=_).*?(?=_)"),
    pca_id = substr(parameter,nchar(parameter),nchar(parameter)),
    predictor = paste0("pc",pca_id),
    # group = interaction(predictor, species, mod, sep = " | "),
    group = interaction(predictor, species, sep = " | "),
    group = factor(group, levels = unique(group[order(predictor, species,decreasing = T)])),
    overlap_zero = lwr <= 0 & upr >= 0)

# Assign species colors
sp.colors_sub <- c("#b5a331","#339d38","#c26a77","#2f2585","#2b695c")

# Create seperate plots for each growth parameter
lapply(c("Linf","g","t"),function(i){
  
  # Subset data
  plot_df <- beta_df %>% filter(response == i) %>% 
    # Create spaces between species
    arrange(desc(predictor), desc(species), desc(mod)) %>%
    mutate(
      group_id = row_number(),
      species_id = as.numeric(group),
      y_pos = group_id + (species_id - 1) * 3
    )
  
  # Plot
  p <- ggplot(data=plot_df, aes(y=y_pos,x=mean,color = species,alpha = overlap_zero)) +
    geom_vline(xintercept = 0, color = "red",linewidth =2) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr),
                   linewidth = 2,height =0)+
    geom_point(aes(shape = mod),
               size = 5)+
    scale_alpha_manual(values = c(`TRUE` = 0.3, `FALSE` = 1),
                       guide = "none") +
    scale_color_manual(values =sp.colors_sub)+
    xlab("")+
    theme(legend.position="none")+
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.border  = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 2),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_text(size = 14)
    )
  print(p)
  
  #File name
  p_name <- paste0("beta_plot_",i,".png")
  
  # Export
  ggsave(
    file.path(export_dir,"beta_tree_plots",p_name),
    plot = p,
    bg = "transparent",
    width = 5,
    height = 10,
    dpi = 300
  )
  
} 
)

# Beta prediction plots (Fig 7)  -----------------------------------------------

# responses used for predictions
pred_vec <- c("ig") 
pred_vec <- c("Linf","ig")   

# X axis limits
xmins <- pred_bridged %>% 
  select(PC1,PC2,PC3) %>% 
  summarise(across(everything(),~min(.x)))
xmaxs <- pred_bridged %>% 
  select(PC1,PC2,PC3) %>% 
  summarise(across(everything(),~max(.x)))

# x axis breaks
x_breaks <- purrr::map2(xmins,xmaxs,function(x,y) round(seq(x,y,length.out =4),2))

# Y axis limits
y_min_vec  <- sapply(pred_vec, function(x) {
  pred_bridged %>% 
    select(contains(x)) %>% 
    summarise(across(everything(),~min(.x))) %>% 
    min()
})
y_max_vec <- sapply(pred_vec, function(x) {
  pred_bridged %>% 
    select(contains(x)) %>% 
    summarise(across(everything(),~max(.x))) %>% 
    max()
})


# Plotting for loop: 
# Each species, response, and predictor receives its own plot
for (j in 1:length(sp)) {
  for (i in 1:length(pred_vec)) {
    for (k in 1:3) {
      
      # Subset data for prediction and CIs
      plot_df <- pred_bridged %>% filter(species == sp[j]) %>% 
        select(contains("PC"),contains(pred_vec[i])) %>% 
        select(contains(as.character(k)))
      colnames(plot_df) <- c("pred","mean","lwr","upr")
      
      
      # Subset data for rug lines
      rug_df <- actual_df %>% 
        filter(species == sp[j]) %>% 
        select(-species,-group_name,-age,-length) %>% 
        distinct() %>% 
        select(contains(as.character(k)))
      colnames(rug_df) <- "rug"
      
      # Set y axis limits based on response
      min_y <- y_min_vec[i]
      max_y <- y_max_vec[i]
      
      # Plotting
      p <- ggplot(data=plot_df,aes(x = pred,y=mean))+
        geom_line(color = sp.colors[j]) +
        geom_ribbon(aes(ymin=lwr, ymax=upr), linewidth=0, alpha=ribbon.alpha,fill = sp.colors[j])+
        geom_rug(data = rug_df,aes(x = rug,size=rug.size),inherit.aes = F,color = sp.colors[j])+
        #ylim(c(min_y,max_y))+
        #xlim(xmins[[k]],xmaxs[[k]])+
        xlab("")+
        ylab("")+
        scale_x_continuous(
          limits = c(xmins[[k]],xmaxs[[k]]),
          breaks = x_breaks[[k]],
          labels = scales::label_number(accuracy = 0.01),
          guide = guide_axis(check.overlap = FALSE))+
        scale_y_continuous(limits = c(min_y,max_y))+
        theme_classic(base_size = base.size)+ 
        theme(legend.position="none") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = border.size)) 
      print(p)
      
      # Create file name based on response, predictor, and species
      beta_plot_name <- paste0(pred_vec[i],"_beta_plots/",sp[j],"_PC",k,"_",pred_vec[i],"_beta_plot.png")
      
      # Export plot
      ggplot2::ggsave(file.path(export_dir,beta_plot_name),
                      p,
                      width = 10,
                      height = 8,
                      dpi = 300)
    }
  }
}

