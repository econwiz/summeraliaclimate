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

#split up hierid into adm0, adm1, adm2... for now. Also make sure to grab all the temp poly
obs <- read_csv(obs_csv, show_col_types = FALSE) %>%
  clean_names() %>%
  separate(hierid, into = c("iso", "adm1_raw", "adm2_raw"),
           sep = "\\.", fill = "right", extra = "drop") %>%
  mutate(
   iso       = toupper(iso),
   adm1_raw  = toupper(adm1_raw),
   adm2_id   = toupper(adm2_raw),
   year      = as.integer(year),
   adm1_id = if_else(
    adm1_raw %in% allowed_ids,
    adm1_raw,
    if_else(
     str_detect(adm1_raw, "^[0-9]+$"),
     str_c(str_sub(iso, 1, 2), str_pad(adm1_raw, 2, pad = "0")),
     str_c(str_sub(iso, 1, 2), adm1_raw)
    )
   )
  )%>%
  select(
    iso, adm1_id, year,
    matches("^tavg_poly_[1-4]_gmfd$", ignore.case = TRUE),
    matches("^lr_tavg_gmfd",          ignore.case = TRUE)
  ) %>%
  distinct(iso, adm1_id, year, .keep_all = TRUE)

#from the in panel, get all the iso, adm1_ids, years. drop any superfluous polynomials from the original panel
panel <- read_dta(panel_in) %>%
  clean_names() %>%
  mutate(
    iso      = toupper(as.character(iso)),
    adm1_id  = toupper(as.character(adm1_id)),
    year     = as.integer(year)
  ) %>%
  select(-matches("gmfd", ignore.case = TRUE), -matches("best", ignore.case = TRUE))

#merge on adm1_id
merged <- left_join(panel, obs, by = c("iso", "adm1_id", "year"))
message("Rows after merge: ", nrow(merged))

gmfd_cols <- grep("gmfd", names(merged), ignore.case = TRUE, value = TRUE)
if (length(gmfd_cols) > 0) {
  names(merged) <- sub("gmfd", "GMFD", names(merged), ignore.case = TRUE)
}

write_dta(merged, panel_out)


