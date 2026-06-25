#-------------------------------------------------------------------------------
#
#  Growth curve model plotting
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: 04-13-2026

# DESCRIPTION: Creates length/growth curve plots, predicted growth at hydrology 
# plots, plots for manuscript appendix.


# House Keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(rstan)
library(parallel)
library(tidyverse)

# Directories
loo_dir <- "loo_outputs_cat"
pred_dir <- "input_data"
plot_dir <- "stan_outputs/plotting_info"
export_dir <- "figures_cat"
mu_dir <- file.path(export_dir,"mu_plots")
fun_dir <- "functions"

# Custom functions
source(file.path(fun_dir,"growth_summary_functions.R"))

# Load data
curve_df <- 
  readRDS(file.path(loo_dir,"stacked_site_curves_2026-06-22.rds"))
cat_curve_df <- 
  readRDS(file.path(loo_dir,"stacked_cat_curves_2026-06-22.rds"))
mu_curve_df <- 
  readRDS(file.path(loo_dir,"stacked_mu_curves_2026-06-22.rds"))
ind_mu_curve_df <- 
  readRDS(file.path(loo_dir,"ind_mu_curves_2026-06-22.rds"))
pred_df <-
  readRDS(file.path(pred_dir,"fsgrw_predictors_2026-06-17.rds"))
age_df <- 
  readRDS(file.path(pred_dir,"fsage_cleaned_2026-06-18.rds"))


# JORFLO DATA  ----------------------------------------------------------------

# Load in jorflo random effect only curve
j_mu_curve_df <- readRDS(
  file.path("loo_outputs","stacked_mu_curves_2026-06-22.rds")
) %>% 
  filter(species == "JORFLO")
j_ind_mu_curve_df <- readRDS(
  file.path("loo_outputs","ind_mu_curves_2026-06-22.rds")
) %>% 
  filter(species == "JORFLO")

# Add jorflo data to curve data.frames
mu_curve_df <- mu_curve_df %>% 
  bind_rows(j_mu_curve_df)
ind_mu_curve_df <- ind_mu_curve_df %>% 
  bind_rows(j_ind_mu_curve_df)


# Prepare data  ----------------------------------------------------------------

# Actual data for points and rug lines
actual_df <- age_df %>% left_join(pred_df) %>% 
  mutate(group_name = paste(region,site,wateryear)) %>% 
  select(species,group_name,age,length,wet_sum_365day,hydroperiod)

# Species names
sp <- unique(actual_df$species)
sp <- sp[order(sp)]

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


# Population curves - model comparisons (Fig s3.1)  ----------------------------

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
  if (sp[i] == "POELAT") {
    mod_compare_sp$top_mod[mod_compare_sp$mod == "lg"] <- 1
    }
  if (sp[i] == "HETFOR") {
    mod_compare_sp$top_mod[mod_compare_sp$mod == "vb"] <- 1
    }
  
  # Assign color scheme
  mod_compare_sp$top_mod <- factor(mod_compare_sp$top_mod,levels = c(0,1))
  top_mod_color <- c("black","red")
  
  # Assign line types for each model
  # Number of models, if 3 opposed to 4, then no stacked model
  n.mods <- n_distinct(mod_compare_sp$mod)
  line_vec <- c(2,3,4,1)[1:n.mods]
  
  # Length at age plot
  length_compare <-ggplot(
    data=mod_compare_sp,
    aes(
      x = age,
      y=length_pred_median,
      group = mod,
      linetype = mod,
      colour = top_mod
      )
    )+
    geom_line() +
    geom_ribbon(
      aes(
        ymin=length_pred_lwr, 
        ymax=length_pred_upr,
        colour = top_mod
        ), 
      linetype=0, 
      alpha=0.1
      )+
    geom_rug(
      data=actual_df %>% filter(species == sp[i]),
      aes(x=age),
      inherit.aes  = F
      )+
    scale_linetype_manual(values =line_vec) +
    scale_color_manual(values =top_mod_color)+
    xlab("")+
    ylab("")+
    xlim(c(0,max_age_mu))+
    theme_classic(base_size = base.size.mu)+ 
    theme(
      legend.position="none",
      panel.border = element_rect(
        color = "black", 
        fill = NA, 
        size = border.size
        )
      )
  print(length_compare)
  
  # Growth at age plot
  growth_compare <- ggplot(
    data=mod_compare_sp,
    aes(
      x = age,
      y=growth_pred_median,
      group = mod,linetype = mod, 
      colour = top_mod
      )
    )+
    geom_line() +
    geom_ribbon(
      aes(
        ymin=growth_pred_lwr, 
        ymax=growth_pred_upr,
        group = mod,
        colour = top_mod
        ), 
      linetype=0, 
      alpha=0.1
      )+
    geom_rug(
      data=actual_df %>% filter(species == sp[i]),
      aes(x=age),
      inherit.aes  = F
      )+
    scale_linetype_manual(values =line_vec) +
    scale_color_manual(values =top_mod_color)+
    xlab("")+
    ylab("")+
    xlim(c(0,max_age_mu))+
    theme_classic(base_size = base.size.mu)+ 
    theme(
      legend.position="none",
      panel.border = element_rect(
        color = "black", 
        fill = NA, 
        size = border.size
        )
      )
  print(growth_compare)
  
  # Export length
  length_compare_name <- paste0("length_age_compare_plot_",sp[i],".png")
  ggplot2::ggsave(
    file.path(mu_dir,length_compare_name),
    length_compare,
    width = 10,
    height = 8,
    dpi = 300
    )
  
  # Export growth
  growth_compare_name <- paste0("growth_age_compare_plot_",sp[i],".png")
  ggplot2::ggsave(
    file.path(mu_dir,growth_compare_name),
    growth_compare,
    width = 10,
    height = 8,
    dpi = 300
    )
  
}



# Population curves - all species (Fig s3.2)  ----------------------------------

# Age at length
mu_l_plot <-ggplot(
  data=mu_curve_df,
  aes(
    x = age,
    y=length_pred_median,
    group = species,
    colour = species)
  )+
  geom_line() +
  geom_ribbon(
    aes(
      ymin=length_pred_lwr, 
      ymax=length_pred_upr,
      fill = species
      ), 
    linetype=0,
    alpha=0.1
    )+
  geom_rug(
    data=actual_df,
    aes(x=age,colour = species),
    inherit.aes  = F)+
  scale_color_manual(values =sp.colors)+
  xlab("")+
  ylab("")+
  xlim(c(0,max_age_mu))+
  theme_classic(base_size = base.size.mu)+ 
  theme(
    legend.position="none"+
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = border.size
      )
    )
ggplot2::ggsave(
  file.path(mu_dir,"mu_length_age_plot.png"),
  mu_l_plot,
  width = 10,
  height = 8,
  dpi = 300
  )


# growth at age
mu_g_plot <-ggplot(
  data=mu_curve_df,
  aes(
    x = age,
    y=growth_pred_median,
    group = species,
    color=species
    )
  )+
  geom_line() +
  geom_ribbon(
    aes(
      ymin=growth_pred_lwr, 
      ymax=growth_pred_upr,
      group = species,
      fill=species
      ), 
    linetype=0, 
    alpha=0.1
    )+
  geom_rug(
    data=actual_df,
    aes(x=age,colour = species),
    inherit.aes  = F
    )+
  scale_color_manual(values =sp.colors)+
  xlab("")+
  ylab("")+
  xlim(c(0,max_age_mu))+
  theme_classic(base_size = base.size.mu)+ 
  theme(
    legend.position="none",
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = border.size)
    )
ggplot2::ggsave(
  file.path(mu_dir,"mu_growth_age_plot.png"),
  mu_g_plot,
  width = 10,
  height = 8,
  dpi = 300
  )



# Hydroperiod specific curves (Fig s3.3) ---------------------------------------
sp_cat <- sp[sp !="JORFLO"]
# Set group colors
cat_col <- c("#132B43","#355E9D","#56B1F7")


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
for(i in 1:length(sp_cat)){
  
  # Subset actual data for ruglines
  l_fish_plot <- actual_df %>% 
    filter(species == sp_cat[i]) 
  
  # Limit plot only predictions within the range of actual fish sampled
  max_age <-plyr::round_any(max(l_fish_plot$age),25,f=ceiling)
  max_age <- 360
  l_plot_df <- cat_curve_df %>% 
    filter(
      species == sp_cat[i],
      age <=max_age
    ) 
  
  # Set y axis limits
  max_length <- plyr::round_any(max(l_plot_df$length_pred_upr),5,ceiling)
  
  # Create plot
  p <-ggplot(
    data=l_plot_df,
    aes(
      x = age,
      y=length_pred_median,
      group = hydroperiod,
      colour = hydroperiod
      )
    )+
    geom_line() +
    geom_ribbon(
      aes(
        ymin=length_pred_lwr, 
        ymax=length_pred_upr,
        fill = hydroperiod
        ), 
      linetype=0, 
      alpha=0.1
      )+
    geom_point(
      data = l_fish_plot,
      aes(
        x=age,
        y = length,
        group = hydroperiod, 
        colour = hydroperiod
        ),
      size = point.size
      )+
    scale_color_manual(values = cat_col)+
    scale_fill_manual(values = cat_col)+
    xlab("")+
    ylab("")+
    scale_y_continuous(
      limits = c(0,max_length),
      #breaks = len_break
    )+
    xlim(c(0,max_age))+
    theme_classic(base_size = base.size)+ 
    theme(
      legend.position="none",
      panel.border = element_rect(
        color = "black", 
        fill = NA, 
        size = border.size
        )
      ) 
  print(p)
  
  # Create file name based on species
  length_plot_name <- paste0(
    "cat_length_at_age_plots/",
    sp_cat[i],
    "_length_age_plot.png"
    )
  
  # Export plot
  ggplot2::ggsave(
    file.path(export_dir,length_plot_name),
    p,
    width = 10,
    height = 8,
    dpi = 300
    )
  
}


# Growth at length plot - species specific
for(i in 1:length(sp_cat)){
  
  # Actual data
  g_fish_plot <- actual_df %>% 
    filter(species == sp_cat[i]) 
  
  # Limit plot only predictions within the range of actual fish sampled
  max_age <-plyr::round_any(max(g_fish_plot$age),25,f=ceiling)
  max_age <- 360
  g_plot_df <- cat_curve_df %>% 
    filter(
      species == sp_cat[i],
      age <= max_age
    )
  
  # Set y axis limits
  max_growth <- plyr::round_any(max(g_plot_df$growth_pred_upr),.1,ceiling)
  max_growth <- max_growth+grow_y_lim_adj[i]
  #growth_break <- seq(0,max_growth,length.out = 4)
  #growth_break <- round(growth_break,2)
  
  # Plotting
  library(ggplot2)
  p <-ggplot(
    data=g_plot_df,
    aes(
      x = age,
      y=growth_pred_median,
      group = hydroperiod,
      color=hydroperiod
      )
    )+
    geom_line() +
    geom_ribbon(
      aes(
        ymin=growth_pred_lwr, 
        ymax=growth_pred_upr,
        group = hydroperiod,
        fill=hydroperiod
        ), 
      linetype=0, 
      alpha=0.1
      )+
    geom_rug(
      data=g_fish_plot,
      aes(x=age,colour = hydroperiod),
      inherit.aes  = F,
      line.width = rug.size
      )+
    scale_color_manual(values = cat_col)+
    scale_fill_manual(values = cat_col)+
    xlab("")+
    ylab("")+
    xlim(c(0,max_age))+
    scale_y_continuous(
      limits = c(0,max_growth),
      #breaks = growth_break
    )+
    theme_classic(base_size = base.size)+ 
    theme(
      legend.position="none",
      panel.border = element_rect(
        color = "black", 
        fill = NA, 
        size = border.size)
      ) 
  print(p)
  
  # Create file name based on species
  growth_plot_name <- paste0(
    "cat_growth_at_age_plots/",
    sp_cat[i],
    "_growth_age_plot.png"
    )
  
  # export
  ggplot2::ggsave(
    file.path(export_dir,growth_plot_name),
    p,
    width = 10,
    height = 8,
    dpi = 300
    )
  
}

# Jorflo plot
j_data <- actual_df %>% filter(species == "JORFLO")
j_length <-ggplot(
  data=j_mu_curve_df,
  aes(x = age,y=length_pred_median)
  )+
  geom_line(color="#97978F") +
  geom_ribbon(
    aes(
      ymin=length_pred_lwr, 
      ymax=length_pred_upr,
      fill = "#97978F"
      ), 
    linetype=0, 
    alpha=0.1
    )+
  geom_point(
    data = j_data,
    aes(x=age,y = length),
    size = point.size,color="#97978F")+
  xlab("")+
  ylab("")+
  scale_y_continuous(
    limits = c(0,max_length),
    #breaks = len_break
  )+
  xlim(c(0,max_age))+
  theme_classic(base_size = base.size)+ 
  theme(
    legend.position="none",
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = border.size
      )
    ) 
print(j_length)

j_growth <-ggplot(
  data=j_mu_curve_df,
  aes(x = age,y=growth_pred_median)
  )+
  geom_line(color = "#97978F") +
  geom_ribbon(
    aes(
      ymin=growth_pred_lwr, 
      ymax=growth_pred_upr,
      fill="#97978F"
      ), 
    linetype=0, 
    alpha=0.1
    )+
  geom_rug(
    data=j_data,aes(x=age),
    inherit.aes  = F,
    line.width = rug.size,
    colour = "#97978F"
    )+
  xlab("")+
  ylab("")+
  xlim(c(0,max_age))+
  scale_y_continuous(
    limits = c(0,max_growth),
    #breaks = growth_break
  )+
  theme_classic(base_size = base.size)+ 
  theme(
    legend.position="none",
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = border.size)
    ) 
print(j_growth)

# Export plot
ggplot2::ggsave(
  file.path(export_dir,"cat_length_at_age_plots/JORFLO_length_age.png"),
  j_length,
  width = 10,
  height = 8,
  dpi = 300
  )
ggplot2::ggsave(
  file.path(export_dir,"cat_growth_at_age_plots/JORFLO_growth_age.png"),
  j_growth,
  width = 10,
  height = 8,
  dpi = 300
  )
