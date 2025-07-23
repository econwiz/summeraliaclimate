#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(haven)
library(janitor)
library(stringr)

#file paths. Panel in is the public mortality dataset merged upon. Obs_csv contains temp/precip
#data. Panel_out is what to write the new merged panel to w/ new obs data 

panel_in   <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
obs_csv    <- "/user/ab5405/summeraliaclimate/code/obs_csvs/ERA5_025_by_region_year.csv"
panel_out  <- "/user/ab5405/summeraliaclimate/code/global_mortality_panel_public_ERA5_025.dta"

#Since the shapefile used already labels countries in the EU, rename all countries in the panel_in
#accordingly, in order to not drop certain observations. 
euro_iso <- c(
  "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP",
  "EST","FIN","GBR","GRC","HRV","HUN","IRL","ISL","ITA",
  "LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
  "POL","PRT","ROU","SVK","SVN","SWE","TUR"
)

#Read the original mortality panel, clean/lowercase the names, and extract iso, adm1_id, and 
#adm_2id. Remove all old temperature data. 
df_panel <- read_dta(panel_in) %>%
  clean_names() %>%
  mutate(
    iso     = str_trim(toupper(iso)),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id = str_trim(toupper(adm1_id)),
    adm2_id = str_trim(toupper(adm2_id)),
    year    = as.integer(year)
  ) %>%
  select(-matches("era5_025|best|gmfd", ignore.case = TRUE))

#Read the new observational dataset, making sure to save 
obs_raw <- read_csv(
  obs_csv,
  col_types = cols(
    iso     = col_character(),
    adm1_id = col_character(),
    adm2_id = col_character(),
    year    = col_integer(),
    .default = col_double()
  ) %>%
  clean_names() %>%
  mutate(
    iso     = str_trim(toupper(iso)),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id = str_trim(toupper(adm1_id)),
    adm2_id = str_trim(toupper(adm2_id)),
    year    = as.integer(year)
  )

#––– 3. merge on ISO + ADM-1 + ADM-2 + year
merged <- df_panel %>%
  left_join(obs_raw, by = c("iso","adm1_id","adm2_id","year"))

#––– 4. uppercase "ERA5_025" in variable names
tidy_names <- function(x) str_replace_all(x, regex("era5-025|era5_025", ignore_case = TRUE), "ERA5_025")
names(merged) <- tidy_names(names(merged))

#––– 5. diagnostics
diag <- merged %>% summarise(
  total_rows    = n(),
  with_climate  = sum(!is.na(tavg_poly_1_ERA5_025)),
  pct_populated = mean(!is.na(tavg_poly_1_ERA5_025)) * 100,
  adm1_clusters = n_distinct(adm1_id),
  adm2_clusters = n_distinct(adm2_id)
)
print(diag)

#––– 6. write out final merged panel
write_dta(merged, panel_out)

