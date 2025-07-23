
#!/usr/bin/env Rscript

library(tidyr)
library(janitor)
library(haven)
library(dplyr)
library(readr)
library(stringr)

# --- Paths ---
panel_in    <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
# point this at your precomputed MERRA-2 by-region/year CSV
obs_csv      <- "/user/ab5405/summeraliaclimate/code/obs_csvs/MERRA2_by_region_year.csv"
panel_out    <- "/user/ab5405/summeraliaclimate/code/regressions/prep_data/global_mortality_panel_public_MERRA2.dta"

# --- 1. Load MERRA-2 observations (already ADM1–year level) ---
merra2_obs <- read_csv(obs_csv, show_col_types = FALSE) %>%
  clean_names() %>%
  # hierid assumed "ISO.ADM1"; adjust regex if different
  extract(hierid, into = c("iso", "adm1_id"), regex = "^([^.]+)\\.([^.]+)", remove = FALSE) %>%
  mutate(
    iso      = toupper(iso),
    adm1_id  = toupper(adm1_id),
    year     = as.integer(year)
  ) %>%
  # drop any hex or malformed ADM1s, if present
  filter(!str_detect(adm1_id, "^0x[0-9A-F]+$")) %>%
  distinct(iso, adm1_id, year, .keep_all = TRUE)

message("✅ Loaded MERRA-2 obs with ", nrow(merra2_obs), " rows (after hex filter)")

# --- 2. Load mortality panel ---
panel <- read_dta(panel_in) %>%
  clean_names() %>%
  mutate(
    iso      = toupper(as.character(iso)),
    adm1_id  = toupper(as.character(adm1_id)),
    year     = as.integer(year)
  ) %>%
  # drop any residual GMFD columns if carried over
  select(-matches("gmfd", ignore.case = TRUE), -matches("best", ignore.case = TRUE))

message("✅ Loaded mortality panel with ", nrow(panel), " rows")

# --- 3. Check for merge mismatches ---
panel_keys <- panel %>% distinct(iso, adm1_id, year)
obs_keys   <- merra2_obs %>% distinct(iso, adm1_id, year)

n_panel    <- nrow(panel_keys)
n_obs      <- nrow(obs_keys)
n_missing  <- anti_join(panel_keys, obs_keys, by = c("iso", "adm1_id", "year")) %>% nrow()

message("🔍 Panel ADM1–year combos:  ", n_panel)
message("🔍 MERRA-2 ADM1–year combos: ", n_obs)
message("⚠️  Panel rows not matched:  ", n_missing)

if (n_missing > 0) {
  warning(n_missing, " panel ADM1–year combos have no MERRA-2 match; they may be dropped in Stata.")
}

# --- 4. Merge ---
merged <- left_join(panel, merra2_obs, by = c("iso", "adm1_id", "year"))
message("📎 Rows after merge: ", nrow(merged))

# --- 5. Rename MERRA-2 variables to uppercase prefix “MERRA2” ---
merra2_cols <- names(merged)[str_detect(names(merged), regex("merra2", ignore_case = TRUE))]
if (length(merra2_cols) > 0) {
  merged <- merged %>%
    rename_with(
      ~ str_replace_all(.x, regex("merra2", ignore_case = TRUE), "MERRA2"),
      .cols = all_of(merra2_cols)
    )
} else {
  message("⚠️  No MERRA-2 columns found in merged panel.")
}

# --- 6. Check for NA in MERRA-2 columns ---
if (length(merra2_cols) > 0) {
  na_rates <- merged %>%
    summarise(across(all_of(str_replace_all(merra2_cols, regex("merra2", ignore_case = TRUE), "MERRA2")),
                     ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing")
  print(na_rates)
}

# --- 7. Save final output ---
write_dta(merged, panel_out)
message("✅ Wrote merged panel to: ", panel_out)

