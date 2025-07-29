#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# This merges MERRA2 observational temperature and precipitation data 
# with the original Carleton et al. (2022) mortality panel.
#
# Input:
# - Original mortality panel (.dta)
# - Aggregated MERRA2-by-region CSV with temperature and precipitation
#
# Output:
# - New panel with MERRA2 vars written to `prep_panels/`
#
# ------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(haven)
library(janitor)
library(stringr)

#get filepaths. Use the original carleton mortality panel to merge on, and select the obs csv for the particular product. Write out a new mortality panel.
panel_in  <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
obs_csv   <- "/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs/MERRA2_by_region_year.csv"
panel_out <- "/user/ab5405/summeraliaclimate/code/regressions/prep_panels/global_mortality_panel_public_MERRA2.dta"

#since this shapefile already codes EU countries as EU, we must recode them in the observational panel first
euro_iso <- c(
  "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP",
  "EST","FIN","GBR","GRC","HRV","HUN","IRL","ISL","ITA",
  "LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
  "POL","PRT","ROU","SVK","SVN","SWE","TUR"
)

#read the original mortality panel, and clean names. Catch any stray EU countries and rename their country code. Drop any old
#weather variables
panel <- read_dta(panel_in) %>%
  clean_names() %>%
  mutate(
    iso      = str_trim(toupper(iso)),
    iso      = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id  = str_trim(toupper(adm1_id)),
    adm2_id  = str_trim(toupper(adm2_id)),
    year     = as.integer(year)
  ) %>%
  select(-matches("merra2|gmfd|best", ignore.case = TRUE))

#read the obs_csv, replace any dashes with underscores, rename EU countries. 
obs_raw <- read_csv(obs_csv, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    iso      = str_trim(toupper(iso)),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id  = str_trim(toupper(adm1_id)),
    adm2_id  = str_trim(toupper(adm2_id)),
    year     = as.integer(year)
  )

# Merge the two panels on their shared codes, and make sure to drop superfluous adm2 lr average from the merge. 
merged <- panel %>%
  left_join(obs_raw, by = c("iso", "adm1_id", "adm2_id", "year"))

#Uppercase all variables
names(merged) <- sub("merra2", "MERRA2", names(merged), ignore.case = TRUE)

#Print diagnostics about the data to see what is missing.
diag <- merged %>% summarise(
  total_rows    = n(),
  with_climate  = sum(!is.na(tavg_poly_1_MERRA2)),
  pct_populated = mean(!is.na(tavg_poly_1_MERRA2)) * 100,
  adm1_clusters = n_distinct(adm1_id),
  adm2_clusters = n_distinct(adm2_id)
)
print(diag)

#save the panels and overwrite original panel if needed.
write_dta(merged, panel_out)

