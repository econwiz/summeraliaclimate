#!/usr/bin/env Rscript

#––– Libraries
library(dplyr)
library(readr)
library(haven)
library(janitor)
library(stringr)

#––– File paths
datapath <- list(
  panel = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta",
  obs   = "/user/ab5405/summeraliaclimate/code/obs_csvs/ERA5_025_by_region_year.csv",
  out   = "/user/ab5405/summeraliaclimate/code/global_mortality_panel_public_ERA5_025.dta"
)

#––– European ISO3 codes to recode to 'EU'
euro_iso <- c(
  "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP",
  "EST","FIN","GBR","GRC","HRV","HUN","IRL","ISL","ITA",
  "LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
  "POL","PRT","ROU","SVK","SVN","SWE","TUR"
)

#––– 1. Read & preprocess mortality panel (keep ADM-2 for FE)
panel <- read_dta(datapath$panel) %>%
  clean_names() %>%
  mutate(
    iso     = str_trim(toupper(iso)),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id = str_trim(toupper(adm1_id)),
    adm2_id = str_trim(toupper(adm2_id)),
    year    = as.integer(year)
  ) %>%
  select(-matches("eras?5[-_]?025|merra2|gmfd|best", ignore.case = TRUE))

#––– 2. Read ERA5-025 observations, convert dashes to underscores in names
raw_obs <- read_csv(datapath$obs, show_col_types = FALSE) 
# replace all dashes in column names with underscores
tidy_names_raw <- names(raw_obs) %>% str_replace_all("-", "_")
names(raw_obs) <- tidy_names_raw

df_obs <- raw_obs %>%
  clean_names() %>%
  mutate(
    iso     = str_trim(toupper(iso)),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id = str_trim(toupper(adm1_id)),
    adm2_id = str_trim(toupper(adm2_id)),
    year    = as.integer(year)
  )

#––– 3. Merge on ISO + ADM-1 + ADM-2 + year
# now using recoded ISO for Europe
df_merged <- panel %>%
  left_join(df_obs, by = c("iso", "adm1_id", "adm2_id", "year"))

#––– 4. Rename ERA5_025 variables to uppercase format
names(df_merged) <- names(df_merged) %>% 
  str_replace_all(regex("era5_025", ignore_case=TRUE), "ERA5_025")

#––– 5. Diagnostics
diag <- df_merged %>% summarise(
  total_rows    = n(),
  with_climate  = sum(!is.na(tavg_poly_1_ERA5_025)),
  pct_populated = mean(!is.na(tavg_poly_1_ERA5_025))*100,
  adm1_clusters = n_distinct(adm1_id),
  adm2_clusters = n_distinct(adm2_id)
)
print(diag)

#––– 6. Write and confirm
write_dta(df_merged, datapath$out)
message("✅ ERA5_025 panel written to: ", datapath$out)

