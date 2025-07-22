#!/usr/bin/env Rscript

# Merge the global mortality panel with GMFD exposure,
# dropping all old GMFD vars from the panel and capitalizing “GMFD”
# in the newly merged columns.

library(haven)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(janitor)

panel_in  <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
obs_csv   <- "/user/ab5405/summeraliaclimate/code/obs_csvs/MERRA2_by_region_year.csv"
panel_out <- "/user/ab5405/summeraliaclimate/code/global_mortality_panel_public_MERRA2.dta"

obs <- read_csv(obs_csv, show_col_types = FALSE) %>%
  clean_names() %>%  # lowercase, hyphens→underscores
  extract(
    hierid,
    into   = c("iso", "adm1_id"),
    regex  = "^([^.]+)\\.([^.]+)",
    remove = FALSE
  ) %>%
  mutate(
    iso      = as.character(iso),
    adm1_id  = as.character(adm1_id),
    year     = as.integer(year)
  ) %>% select(
    iso, adm1_id, year,
    starts_with("tavg_poly"),        # temp polys
    starts_with("lr_tavg"),          # linear-response term(s)   # <— NEW: precip polynomials
  ) %>%
  distinct(iso, adm1_id, year, .keep_all = TRUE) 

message("MERRA2 columns in obs:\n  ", 
        paste(names(obs)[str_detect(names(obs), "merra2")], collapse = "\n  "))

panel <- read_dta(panel_in) %>%
  mutate(
    iso      = as.character(iso),
    adm1_id  = as.character(adm1_id),
    year     = as.integer(year)
  )%>% select(-matches("gmfd", ignore.case = TRUE))

panel_keys <- panel %>% distinct(iso, adm1_id, year)
obs_keys   <- obs   %>% distinct(iso, adm1_id, year)

n_panel   <- nrow(panel_keys)
n_obs     <- nrow(obs_keys)
n_missing <- panel_keys %>%
  anti_join(obs_keys, by = c("iso", "adm1_id", "year")) %>%
  nrow()

message("Distinct ADM1–year in panel: ", n_panel)
message("Distinct ADM1–year in obs:   ", n_obs)
message("Panel keys not in obs:      ", n_missing)
if (n_missing > 0) {
  warning(n_missing,
          " ADM1–year keys lack GMFD data; Stata prep will drop them.")
}

merged <- panel %>%
  left_join(obs, by = c("iso", "adm1_id", "year"))

message("Rows before merge: ", nrow(panel))
message("Rows after  merge: ", nrow(merged))

merged <- merged %>%
  rename_with(
    .cols = matches("merra2", ignore.case = TRUE),
    .fn   = ~ str_replace_all(.x, regex("merra2", ignore_case = TRUE), "MERRA2")
  )

m2_cols <- names(merged)[str_detect(names(merged), "MERRA2")]
na_rates <- merged %>%
  summarise(across(all_of(m2_cols), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing")
print(na_rates)

write_dta(merged, panel_out)
message("✅ Wrote merged panel to: ", panel_out)

