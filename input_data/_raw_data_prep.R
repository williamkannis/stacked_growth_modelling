#-------------------------------------------------------------------------------
#
#  Raw data preparation          
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: Feb 2, 2026

# DESCRIPTION: Prepares raw data to a state that can be shared with manuscript
# and used to replicate results. This involves coarsening the throw trap and
# hydrology data to annual means. Full age-at-length is shared, and here column
# names are cleaned up and unnecassary fields are removed. No data imputation
# or filtering takes place at this stage.


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages  
library(dplyr)
library(readxl)
library(janitor)
library(purrr)

# Directories
age_dir <- "_raw_data"
raw_dir <- paste0(
  "~/Documents/Work/Everglades post-doc/",
  "Data analysis/Data cleaning/cleaned_data"
)
input_dir <-"input_data"

# Data
age_df <- read_excel(file.path(age_dir,"gatto_age_length.xlsx"))
pis_df <- readRDS(file.path(raw_dir,"pisc_cleaned_2026-02-25.rds"))
len_df <- readRDS(file.path(raw_dir,"fslen_cleaned_2026-02-25.rds"))
phy_df <- readRDS(file.path(raw_dir,"phys_cleaned_2026-02-25.rds"))


# Clean age data ---------------------------------------------------------------

# Data contains some lengths for males, but this is not consistent across 
#  sampling events or specie and can bias growth estimates.

age_clean <- age_df %>% 
  
  # Clean up names
  clean_names() %>% 
  rename(age = ring_count) %>% 
  
  # Fix site names to match main data
  mutate(
    site = case_when(
      site == "07" & region == "SRS" ~ "7",
      T ~ site
    ),
  ) %>% 
  
  # Select relevant fields
  select(wateryear,region,site,species,sex,age,length,weight)


# Final check
summary(age_clean)


# Species key and mean lengths for density data  -------------------------------

# Mean length used for growth at typical length predictions
mean_len_df <- len_df %>% 
  filter(species %in% unique(age_clean$species)) %>% 
  group_by(species) %>% 
  summarise(length = mean(length,na.rm = T)) %>% 
  right_join(key) %>% 
  select(species,sci_name,sci_name_abv,length)


# Growth sampling periods  -----------------------------------------------------

# Create data frame with the wateryear, period, region, site (i.e., sampling
# period) of each growth estimate

samp_df <- age_clean %>% 
  mutate(period = 4) %>% 
  distinct(wateryear,period,region,site) 

# age data were collected during October, which means fish were collected in
# the middle of a water year. To have an accurate one year lag in fish/physical
# predictors we need to create a new year column for the growth measures. Hydro
# data has annual lags for each period, so this is okay. Some data were 
# collected in calender year 2025 but this is not enough for a growth year so 
# remove 2025 from all annual estimations for now

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


# Some sampling events in the age data did not have elctrofishing conducted at
# their respective sample site.For these, impute pisc_index of using the median
# value for that region and year in the published scripts.


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

setdiff(
  fsden_per %>% distinct(wateryear,region,site),
  fsden_year %>% distinct(wateryear,region,site)
) %>% 
  print(n=50)
setdiff(
  fsden_year %>% distinct(wateryear,region,site),
  fsden_per %>% distinct(wateryear,region,site)
) %>% print(n=50)
# The period data doesn't have water 2025 and annual doesn't have 1995, this is 
# okay as we it is an artifact of how the data were collected starting the 
# start of the 1996 calender year, half between the wateryear. we are not using 
# these data fro growth so this is okay, but will need to think about this in 
# the future


# Hydrology  -------------------------------------------------------------------

# Summarize annual hydrological data at the site year for each water period to 
# estimate the impacts of hydrological disturbance and energy on fish growth.
# Also creates categorical groupings based on annual hydroperiod (i.e. days 
# flooded).

hydro_df <- phy_df %>% 
  filter(period == 4) %>% 
  group_by(wateryear,period,region,site) %>% 
  summarise(
    depth_ave_365day = mean(depth_ave_365day,na.rm=T),
    wet_sum_365day = mean(wet_sum_365day,na.rm=T),
    dsldd = mean(dsldd,na.rm=T)
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
comp_df <- hydro_df %>% 
  left_join(pis_df_for) %>% 
  left_join(fsden_year) %>% 
  
  # add column denoting samples iwth age data available
  left_join(samp_df %>% mutate(age_data = T))
  
# Check for collinearity in data
comp_df %>% 
  select(-wateryear,-period,-region,-site,-hydroperiod,-age_data) %>% 
  cor(use="complete.obs")

# Hydrological and pisc data are highly correlated. PCA should be conducted to
# reduce dimensionality


# Export data.frame  -----------------------------------------------------------

# Age-at-length
saveRDS(age_clean,file.path(input_dir,"FCE1302_fsage_at_length.rds"))

# Mean length data
write.csv(
  mean_len_df,
  file.path(input_dir,"FCE1302_fskey_meanlen.csv"),
  row.names = F
  )

# Predictor data
saveRDS(comp_df,file.path(input_dir,"FCE1302_fsgrw_predictors.rds"))

