#-------------------------------------------------------------------------------
#
#  Growth predictor data preparation 
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: Feb 26, 2026

# DESCRIPTION: Compiles and formats data for use as predictor variables for growth
# parameters in hierarchical growth models. The Predictor variables were highly 
# inter-correlated, so we performed principle components analysis to reduce 
# predictors into three composite variables. We also provide hydroperiod (i.e, #
# of days flooded) classifications for use in categorical second-level effects.


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages
library(dplyr)
library(vegan)

# Directories
data_dir <- "~/Documents/Work/Everglades post-doc/Data analysis/Data cleaning/cleaned_data"
input_dir <- "input_data"

# Data
age_df <- readRDS(file.path(input_dir,"fsage_cleaned_2026-02-26.rds"))
pis_df <- readRDS(file.path(data_dir,"pisc_cleaned_2026-02-25.rds"))
len_df <- readRDS(file.path(data_dir,"fslen_cleaned_2026-02-25.rds"))
phy_df <- readRDS(file.path(data_dir,"phys_cleaned_2026-02-25.rds"))


# Growth sampling periods  -----------------------------------------------------

# Create data frame with the wateryear, period, region, site (i.e., sampling
# period) of each growth estimate

samp_df <- age_df %>% 
  distinct(wateryear,period,region,site) 

# age data were collected during October, which means fish were collected in
# the middle of a water year. To have an accurate one year lag in fish/physical
# predictors we need to create a new year column for the growth measures. Hydro
# data has annual lags for each period, so this is okay. Some data were collected
# in calender year 2025 but this is not enough for a growth year so remove 2025
# from all annual estimations for now

grow_year <- phy_df %>% 
  distinct(wateryear,year,period,region,site) %>% 
  mutate(
    growth_year = case_when(
      period %in% 1:4 ~ year,
      period == 5 ~ year+1
    )
  ) %>% 
  select(-year) 

# does growth_year contain all sites in length and physical data?
len_df %>% 
  anti_join(grow_year,join_by(wateryear,period,region,site)) %>% 
  nrow()==0
phy_df %>% 
  anti_join(grow_year,join_by(wateryear,period,region,site)) %>% 
  nrow()==0


# Piscivore CPUE  --------------------------------------------------------------

# Piscivore catch per unit effort data to measure the effects of predator 
# presence on fish growth. Year column needs to be renamed water year column to
# merge in with growth data. This is okay because electrofishing only occurred 
# during the wet season, so year of index will correspond to water year

# TSL data are coarsened to site level, but does not include the sh divsions.

# duplicate the TSL MD and TS into sh sites in hydro data
name_ids <-phy_df %>% 
  distinct(wateryear,region,site) %>% 
  rename(site_full = site) %>% 
  mutate(
    site = case_when(
      substr(site_full,1,2) == "MD" ~ "MD",
      substr(site_full,1,2) == "TS" ~ "TS",
      T~site_full
    )
  )
pis_df_for <- pis_df %>% 
  rename(wateryear = year) %>% 
  filter(region != "PHD") %>% 
  full_join(name_ids) %>% 
  select(-site) %>% 
  rename(site = site_full)
  
  
# Some sampling events in the age data did not have elctrofishing conducted,
# likely due to low water conditions. For these, impute pisc_index of zero.
# Many missing values exist outside of the age sampling events and these
# will not be imputed for pca data. This decision does not change clusters or 
# PCA axes much.

# Impute missing pisc data from age sampling events
pis_df_age <- pis_df_for %>% 
  right_join(samp_df) %>% 
  mutate(
    pisc_index = case_when(
      is.na(pisc_index) ~ 0,
      T~pisc_index
    )
  ) %>% 
  select(-period)

# Extract pisc data that is not linked to age sampling events
pis_df_sub <-
  pis_df_for %>% 
  anti_join(samp_df) %>% 
  filter(!is.na(pisc_index))

# Merge in imputed pisc data with the non imputed pisc data
pis_df_all <- pis_df_age %>% 
  bind_rows(pis_df_sub) %>% 
  filter(!region %in% c("PHD"))


# Fish density  ----------------------------------------------------------------

# estimate the density of fish at each site, each period and each year as a 
# measure of the effects of competition on fish growth

# Period level density
fsden_per <- len_df %>% 
  group_by(wateryear,region,site,period) %>% 
  summarise(
    n_fish = length(species[species != "NOFISH"]),
    area = n_distinct(plot,throw)
    ) %>% 
  ungroup() %>% 
  mutate(fsden_period = n_fish/area) %>% 
  select(-n_fish,-area)

# Annual mean density
fsden_year <-fsden_per %>% 
  left_join(grow_year) %>% 
  group_by(growth_year,region,site) %>% 
  summarise(fsden_annual = mean(fsden_period,na.rm=T))%>% 
  ungroup() %>% 
  rename(wateryear=growth_year) %>% 
  filter(wateryear != 2025)  # filter out 2025 for now. no october data yer

# Does annual mean density have the same rows as year/sites in period level 
# fishdensity?
nrow(fsden_year) == fsden_per %>% distinct(wateryear,region,site) %>% nrow()

fsden_year %>% 
  group_by(wateryear,region,site) %>% 
  filter(n()>1)

setdiff(fsden_per %>% distinct(wateryear,region,site),fsden_year %>% distinct(wateryear,region,site)) %>% print(n=50)
setdiff(fsden_year %>% distinct(wateryear,region,site),fsden_per %>% distinct(wateryear,region,site)) %>% print(n=50)
# The period data doesnt hav water 2025 and naanual doesnt have 1995, this is okay as we it is an artifact of
# how the data were collected starting the start of the 1996 calender year, half betweent he wateryear.
# we are not using these data fro growth so this is okay, but will need to think about
# this in the future

# Hydrology  -------------------------------------------------------------------

# Summarize annual hydrological data at the site year for each water period to 
# estimate the impacts of hydrological disturbance and energy on fish growth.
# Also creates categorical groupings based on anual hydroperiod (i.e. days 
# flooded).

hydro_df <- phy_df %>% 
  group_by(wateryear,period,region,site) %>% 
  summarise(
    depth_ave_365day = mean(depth_ave_365day,na.rm=T),
    wet_sum_365day = mean(wet_sum_365day,na.rm=T),
    dsldd = mean(dsldd,na.rm=T),
    lastdaydry = mean(lastdaydry,na.rm=T)
  ) %>% 
  ungroup() %>% 
  mutate(
    hydroperiod = case_when(
      wet_sum_365day > 360 ~ "long",
      wet_sum_365day >=320 & wet_sum_365day <= 360 ~ "intermediate",
      wet_sum_365day < 320 ~ "short",
      T~NA
    ),
    hydroperiod = factor(hydroperiod,levels=c("short","intermediate","long"))
  )


# Compile all predictors into one data frame  ----------------------------------
comp_df <- age_df %>% 
  distinct(wateryear,period,region,site) %>% 
  left_join(pis_df_all) %>% 
  left_join(fsden_per) %>% 
  left_join(fsden_year) %>% 
  left_join(hydro_df) 

# Check for collinearity in data
comp_df %>% 
  select(-wateryear,-period,-region,-site,-hydroperiod) %>% 
  cor(use="complete.obs")


# Hydrological and pisc data are highly correlated. PCA should be conducted to
# reduce dimensionality


# Create composite variables with PCA  -----------------------------------------

# Use full pisc and hydrology datas ets to conduct PCA on biological and
# hydrology variables

# Prepare data input
pca_data <-hydro_df %>% 
  filter(period == 4) %>% # only use data from period 4 (i.e., when age is collected)
  right_join(pis_df_all) %>% # only use hydrology data that isn't missing pisc data
  left_join(fsden_year)  # load in annual fish densities
pca_input <- pca_data %>% select(-wateryear,-period,-region,-site,-hydroperiod)
  
# Run pca
pca_out <- rda(pca_input,scale = T)

# Examine results
summary(pca_out)
pca_out$CA$v

# Extract pcs that explain atleast 75% of varaition and create dataframe with
# sample event identifiers
n_axes <- 3
pca_id <- hydro_df %>% select(wateryear,period,region,site) %>% filter(period==4) %>% 
  right_join(pis_df_all) %>%select(-pisc_index)
# pca_result <- pca_out$CA$u[,1:n_axes] 
pca_result <-scores(pca_out,choices = c(1,2,3),display = "sites")
pca_df <- cbind(pca_id,pca_result)


# Create and export final predictor output data frame  -------------------------
pred_df <- comp_df %>% 
  left_join(pca_df)

# Export
saveRDS(pred_df,file.path(input_dir,paste0("fsgrw_predictors_",Sys.Date(),".rds")))


# PCA plotting (Fig 2) ---------------------------------------------------------

## PCA plot ##

# Extract scores
site_score_df <- scores(pca_out,choices = c(1,2,3),display = "sites")

# Add identify for sites in growth study
plot_id <- pca_id %>% 
  left_join(comp_df %>% select(wateryear,region,period,site) %>% mutate(in_study = 1)) %>% 
  mutate(in_study = case_when(
    is.na(in_study) ~ 0,
    T~ in_study
  ),
  in_study = factor(in_study))
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
  pca_file <- paste0("pca_plot_",pca_ax[1,i],"-",pca_ax[2,i],".png")
  
  # Plot
  p <-ggplot(data = point_df, aes(x= x,y=y,colour = in_study))+
    geom_point(size=3)+
    scale_color_manual(values = c("grey","black"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_text(data = label_df,aes(x=x,y=y,label=label),inherit.aes = F,color = "darkblue")+
    xlim(xlims[[i]])+
    xlab("")+
    ylab("")+
    theme_classic(base_size = 30)+
    theme(panel.border = element_rect(color = "black", fill = NA, size = 3)) +
    theme(legend.position="none")
  print(p)
  
  ggsave(filename = file.path("figures/pca_plot",pca_file),
         plot = p,
         width = 10,
         height = 8,
         dpi = 300)
}

## Correlation plot ##

# Create cor matrix
pca_cor <- pca_df %>% 
  left_join(pca_data) %>% 
  select(-wateryear,-period,-region,-site,-hydroperiod) %>% 
  cor(use="complete.obs") %>% 
  as.matrix()
pca_cor <- pca_cor[-c(1:3),1:3]

# Export plot
corrplot::corrplot(pca_cor)

png("figures/corr_plot.png", width = 900, height = 1800)
corrplot::corrplot(
  corr = pca_cor,
  cl.pos = "n",
  tl.pos = "n"
)
dev.off()
