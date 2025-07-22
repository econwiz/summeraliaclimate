#!/usr/bin/env Rscript

library(tidyr)
library(janitor)
library(haven)
library(dplyr)
library(readr)
library(stringr)

panel_in   <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
obs_csv    <- "/user/ab5405/summeraliaclimate/code/obs_csvs/GMFD_by_region_year.csv"
panel_out  <- "/user/ab5405/summeraliaclimate/code/global_mortality_panel_public_GMFD.dta"

#pull out the unique adm1 ids from the in panel to match new obs data
allowed_ids <- read_dta(panel_in) %>%
  mutate(adm1_id = toupper(as.character(adm1_id))) %>%
  pull(adm1_id) %>%
  unique()

#split up hierid into adm0, adm1, adm2... skip hex values? for now. Also make sure to grab all the temp polys
obs <- read_csv(obs_csv, show_col_types = FALSE) %>%
  clean_names() %>%
  separate(hierid, into = c("iso", "adm1_raw", "adm2_raw"), sep = "\\.", fill = "right") %>%
  mutate(
    iso       = toupper(iso),
    adm1_raw  = toupper(adm1_raw),
    adm2_raw  = toupper(adm2_raw),
    year      = as.integer(year)
  ) %>%
  select(
    iso, adm1_id, year,
    matches("^tavg_poly_[1-4]_gmfd$", ignore.case = TRUE),
    matches("^lr_tavg_gmfd",        ignore.case = TRUE)
  ) %>% 
  distinct(iso, adm1_id, year, .keep_all = TRUE)

names(obs)
#from the in panel, get all the iso, adm1_ids, years. drop any superfluous polynomials from the original panel
panel <- read_dta(panel_in) %>%
  clean_names() %>%
  mutate(
    iso      = as.character(iso),
    adm1_id  = as.character(adm1_id),
    year     = as.integer(year)
  ) %>%
  select(
    everything(),
    -matches("gmfd", ignore.case = TRUE),
    -matches("best", ignore.case = TRUE)
  )

merged <- panel %>%
  left_join(obs, by = c("iso", "adm1_id", "year"))
message("Original panel rows: ", nrow(panel))
message("Merged panel rows:   ", nrow(merged))

new_cols <- grep("gmfd", names(merged), ignore.case = TRUE, value = TRUE)
if (length(new_cols) > 0) {
  merged <- merged %>%
    rename_with(
      .cols = all_of(new_cols),
      .fn   = ~ str_replace_all(.x, regex("gmfd", ignore_case = TRUE), "GMFD")
    )
}

if (length(new_cols) > 0) {
  gmfd_cols_after <- grep("GMFD", names(merged), value = TRUE)
  na_rates <- merged %>%
    summarise(across(all_of(gmfd_cols_after), ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing")
  print(na_rates)
}


write_dta(merged, panel_out)
message("✅ Merged panel written to: ", panel_out)

