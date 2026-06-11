#-------------------------------------------------------------------------------
#
#  Age at length data cleaning          
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis
# CREATED: Feb 2, 2026

# DESCRIPTION: Combines otilith data with annual hydroperiod for each site. 
# Combined data is then cleaned, removing any male's and missing data.


# Housekeeping  ----------------------------------------------------------------
rm(list = ls())

# Load in packages  
library(dplyr)
library(readxl)
library(janitor)

# Directories
phy_dir <- "~/Documents/Work/Everglades post-doc/Data analysis/Data cleaning/cleaned_data"
raw_dir <- "raw_data"
input_dir <-"input_data"

# Data
age_length_df <- read_excel(file.path(raw_dir,"gatto_age_length.xlsx"))
phy_df <- readRDS(file.path(phy_dir,"phys_cleaned_2026-02-25.rds"))


# Clean data -------------------------------------------------------------------

# Data contains some lengths for males, but this is not consistent across sampling
# events or specie and can bias growth estimates. Here, males are removed. In 
# in addition other missing data are removed

age_length_clean <- age_length_df %>% 
  clean_names() %>% 
  
  # Fix site names to match main data
  mutate(
    site = case_when(
      site == "07" & region == "SRS" ~ "7",
      T ~ site
    ),
    # add period to merge with hydro data
    period = 4  # samples collected in October
  ) %>% 
  
  # Remove males and missing data
  filter(sex != "M",
       !is.na(ring_count),
       !is.na(length))  %>% 
  select(wateryear,period,region,site,species,sex,ring_count,length)


# Final check  -----------------------------------------------------------------
summary(age_length_clean)

# Are all males removed
age_length_clean %>% filter(sex == "M") %>% nrow() ==0

# NAs?
age_length_clean %>% filter(is.na(ring_count)) %>% nrow() ==0
age_length_clean %>% filter(is.na(length)) %>% nrow() ==0


# Export  ----------------------------------------------------------------------
saveRDS(age_length_clean,
        file.path(input_dir,paste0("fsage_cleaned_",Sys.Date(),".rds")))

