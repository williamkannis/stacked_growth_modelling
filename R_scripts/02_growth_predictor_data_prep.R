#-------------------------------------------------------------------------------
#
#  Growth analysis data preparation 
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: Feb 26, 2026

# DESCRIPTION: Prepares data for growth analysis. Cleans otolith-derived 
# age-at-length data, removing any male's and missing data. Data for growth 
# predictors were highly inter-correlated, so we performed principle components 
# analysis to  reduce predictors into three composite variables


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)
library(vegan)
library(ggplot2)

# Directories
input_dir <- "input_data"
fig_dir <- "figures"

# Data
age_df <- readRDS(file.path(input_dir,"FCE1302_fsage_at_length.rds"))
pred_df <- readRDS(file.path(input_dir,"FCE1302_fsgrw_predictors.rds"))


# Age data prep  ---------------------------------------------------------------

# Combine age and predictor data.frames
age_filter <- age_df %>% 
  
  # Remove male and NA lengths from age data
  filter(
    sex != "M",
    !is.na(age),
    !is.na(length)
  )

# Export for analysis
saveRDS(age_filter,file.path(input_dir,"fsage_filtered.rds"))


# Growth sampling periods  -----------------------------------------------------

# Create data frame with the wateryear, period, region, site (i.e., sampling
# period) of each growth estimate

samp_df <- age_filter %>% 
  distinct(wateryear,region,site) 


# Piscivore CPUE data imputation  ----------------------------------------------


# Some sampling events in the age data did not have elctrofishing conducted at
# their respective sample site.For these, impute pisc_index of using the median
# value for that region and year.
# Many missing values exist outside of the age sampling events and these
# will not be imputed for pca data. This decision does not change clusters or 
# PCA axes much.

# Impute missing pisc data from age sampling events using average values for
# that region and year
pisc_avg <- pred_df %>% 
  group_by(region,wateryear) %>% 
  summarise(pisc_index_impute = mean(pisc_index,na.rm=T))

pred_impute <- pred_df %>% 
  left_join(pisc_avg) %>% 
  mutate(
    pisc_index = case_when(
      is.na(pisc_index) & age_data == T ~ pisc_index_impute,
      T~pisc_index
    )
  ) %>% 
  select(-pisc_index_impute)


# Check for collinearity in data -----------------------------------------------
pred_impute %>% 
  select(-wateryear,-period,-region,-site,-hydroperiod,-age_data) %>% 
  cor(use="complete.obs")


# Hydrological and pisc data are highly correlated. PCA should be conducted to
# reduce dimensionality


# Create composite variables with PCA  -----------------------------------------

# Use full pisc and hydrology datas ets to conduct PCA on biological and
# hydrology variables

# Prepare data input
pca_data <-pred_impute %>% 
  filter(!is.na(pisc_index))
pca_input <- pca_data %>% 
  select(
    depth_ave_365day,
    wet_sum_365day,
    dsldd,
    pisc_index,
    fsden_annual)
  
# Run pca
pca_out <- rda(pca_input,scale = T)

# Examine results
summary(pca_out)
pca_out$CA$v

# Extract pcs that explain atleast 75% of variation and create data.frame with
# sample event identifiers
n_axes <- 3
pca_id <- pca_data %>% 
  select(wateryear,region,site) 
pca_result <-scores(pca_out,choices = c(1,2,3),display = "sites")
pca_df <- cbind(pca_id,pca_result)


# Create and export final pca output data frame  -------------------------------
out_df <- samp_df %>% 
  left_join(pca_df)

# Export
saveRDS(out_df,file.path(input_dir,"fsgrw_pca_out.rds"))


# PCA plotting (Fig 2) ---------------------------------------------------------

## PCA plot ##

# Extract scores
site_score_df <- scores(pca_out,choices = c(1,2,3),display = "sites")

# Add identify for sites in growth study
plot_id <- pca_id %>% 
  left_join(pred_df) %>% 
  mutate(age_data = case_when(
    is.na(age_data) ~ 0,
    T~ age_data
  ),
  in_study = factor(age_data))
pca_plot_df <- cbind(plot_id,site_score_df)


# Format species labels
sp_scores_df <- scores(pca_out,choices = c(1,2,3),display = "species") %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("label") %>% 
  mutate(label = case_when(
    label == "fsden_annual" ~ "Fish density",
    label == "pisc_index" ~ "Predator CPUE",
    label == "dsldd" ~ "DSD",
    label == "lastdaydry" ~ "Dry Length",
    label == "wet_sum_365day" ~ "Hydroperiod",
    label == "depth_ave_365day" ~ "Depth"
  ))

# Plot parameters
pca_ax <- combn(c("PC1","PC2","PC3"),2)
xlims <- list(c(-2.5,3),c(-2.5,3),c(-2.5,3))

# Plotting loop
for (i in 1:ncol(pca_ax)) {
  
  # Load data
  point_df <- pca_plot_df
  label_df <- sp_scores_df
  
  # Create x and y axes
  point_df[,"x"] <- point_df[pca_ax[1,i]]
  point_df[,"y"] <- point_df[pca_ax[2,i]]
  label_df[,"x"] <- label_df[pca_ax[1,i]]
  label_df[,"y"] <- label_df[pca_ax[2,i]]
  
  # file name
  pca_file <- paste0("_pca_plot_",pca_ax[1,i],"-",pca_ax[2,i],".png")
  
  # Plot
  p <-ggplot(
    data = point_df, 
    aes(x= x,y=y,colour = in_study)
    )+
    geom_point(size=3)+
    scale_color_manual(values = c("grey","black"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_text(
      data = label_df,
      aes(x=x,y=y,label=label),
      inherit.aes = F,
      color = "darkblue"
      )+
    xlim(xlims[[i]])+
    xlab("")+
    ylab("")+
    theme_classic(base_size = 30)+
    theme(
      legend.position="none",
      panel.border = element_rect(
        color = "black", 
        fill = NA, 
        size = 3
        )
      )
  print(p)
  
  ggsave(
    filename = file.path(
      fig_dir,
      "figure_2",
      pca_file),
     plot = p,
     width = 10,
     height = 8,
     dpi = 300)
}

## Correlation plot ##

# Create cor matrix
pca_cor <- pca_df %>% 
  left_join(pca_data) %>% 
  select(-wateryear,-period,-region,-site,-hydroperiod,-age_data) %>% 
  cor(use="complete.obs") %>% 
  as.matrix()
pca_cor <- pca_cor[-c(1:3),1:3]

# Export plot
corrplot::corrplot(pca_cor)

png(
  file.path(
    fig_dir,
    "figure_2",
    "_corr_plot.png"
    ), 
  width = 900, 
  height = 1800
  )
corrplot::corrplot(
  corr = pca_cor,
  cl.pos = "n",
  tl.pos = "n"
)
dev.off()

