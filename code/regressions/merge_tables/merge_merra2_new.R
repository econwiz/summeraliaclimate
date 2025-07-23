#!/usr/bin/env Rscript

#––– Libraries
library(dplyr)
library(readr)
library(haven)
library(janitor)
library(stringr)

#––– 1. file paths
panel_in  <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
obs_csv   <- "/user/ab5405/summeraliaclimate/code/regressions/generate_obs_data/obs_csvs/MERRA2_by_region_year.csv"
panel_out <- "/user/ab5405/summeraliaclimate/code/regressions/prep_data/global_mortality_panel_public_MERRA2.dta"

#––– 2. European ISOs to recode to 'EU'
euro_iso <- c(
  "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP",
  "EST","FIN","GBR","GRC","HRV","HUN","IRL","ISL","ITA",
  "LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
  "POL","PRT","ROU","SVK","SVN","SWE","TUR"
)

#––– 3. read & prep mortality panel (keep ADM-2 for FE)
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

#––– 4. read & prep GMFD observations
obs_raw <- read_csv(obs_csv, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    iso      = str_trim(toupper(iso)),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id  = str_trim(toupper(adm1_id)),
    adm2_id  = str_trim(toupper(adm2_id)),
    year     = as.integer(year)
  )

#––– 5. merge on ISO + ADM-1 + ADM-2 + year (full ADM-2 alignment)
merged <- panel %>%
  left_join(obs_raw, by = c("iso", "adm1_id", "adm2_id", "year"))

#––– 6. tidy variable names (GMFD → uppercase)
names(merged) <- sub("merra2", "MERRA2", names(merged), ignore.case = TRUE)

#––– 7. diagnostics
diag <- merged %>% summarise(
  total_rows    = n(),
  with_climate  = sum(!is.na(tavg_poly_1_MERRA2)),
  pct_populated = mean(!is.na(tavg_poly_1_MERRA2)) * 100,
  adm1_clusters = n_distinct(adm1_id),
  adm2_clusters = n_distinct(adm2_id)
)
print(diag)

#––– 8. write out final merged panel
write_dta(merged, panel_out)

