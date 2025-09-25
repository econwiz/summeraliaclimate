#!/usr/bin/env Rscript
# Merge ANY *_by_region_year.csv into the Carleton mortality panel
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(haven)
  library(janitor); library(stringr); library(tools)
})

# --------- USER SETTINGS (edit paths if needed) -----------------
panel_in <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
obs_dir  <- "/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs"
out_dir  <- "/user/ab5405/summeraliaclimate/code/regressions/prep_panels"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# If you want to restrict to specific products, put them here (otherwise it auto-detects):
# examples: c("ERA5_025","GMFD","JRA_3Q","MERRA2")
only_products <- c("MERRA2")
# ----------------------------------------------------------------

# EU collapse set (match Kevin’s convention incl. FRA, LIE, TUR)
euro_iso <- c(
  "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP",
  "EST","FRA","FIN","GBR","GRC","HRV","HUN","IRL","ISL","ITA",
  "LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
  "POL","PRT","ROU","SVK","SVN","SWE","TUR"
)

# ---- Load & normalize panel (drop any old weather vars) ----
panel <- read_dta(panel_in) %>%
  clean_names() %>%
  mutate(
    iso     = toupper(str_trim(as.character(iso))),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id = str_trim(as.character(adm1_id)),
    adm2_id = str_trim(as.character(adm2_id)),
    year    = as.integer(year)
  ) %>%
  select(-matches("era5[-_]?025|merra2|gmfd|jra|best", ignore.case = TRUE))

# ---- Find obs CSVs ----
obs_files <- list.files(obs_dir, pattern = "_by_region_year\\.csv$", full.names = TRUE)

if (length(only_products) > 0) {
  # Keep only files whose inferred product matches requested set
  keep <- vapply(obs_files, function(p) {
    # infer product like "ERA5_025" from filename
    base <- basename(p)
    prod <- toupper(gsub("-", "_", sub("_by_region_year\\.csv$", "", base)))
    prod %in% toupper(gsub("-", "_", only_products))
  }, logical(1))
  obs_files <- obs_files[keep]
}

stopifnot(length(obs_files) > 0)

# Helper to infer product tag from filename, e.g. ERA5_025
infer_product <- function(path) {
  nm <- basename(path)
  toupper(gsub("-", "_", sub("_by_region_year\\.csv$", "", nm)))
}

# Keys
key <- c("iso","adm1_id","adm2_id","year")

for (obs_csv in obs_files) {
  product_tag <- infer_product(obs_csv)   # e.g., "ERA5_025", "GMFD", "JRA_3Q", "MERRA2"
  message("\n==== Merging product: ", product_tag, " ====")
  message("Obs file: ", obs_csv)

  # Read & normalize obs
  obs_raw <- read_csv(obs_csv, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(
      iso     = toupper(str_trim(as.character(iso))),
      iso     = if_else(iso %in% euro_iso, "EU", iso),
      adm1_id = str_trim(as.character(adm1_id)),
      adm2_id = str_trim(as.character(adm2_id)),
      year    = as.integer(year)
    )

  # --- Key diagnostics BEFORE merge ---
  panel_keys <- panel %>% distinct(across(all_of(key)))
  obs_keys   <- obs_raw %>% distinct(across(all_of(key)))

  p_not_o <- anti_join(panel_keys, obs_keys, by = key)
  o_not_p <- anti_join(obs_keys, panel_keys, by = key)

  cat("\n=== KEY CHECK (before merge) ===\n")
  cat("Panel keys: ", nrow(panel_keys), "\n")
  cat("Obs keys:   ", nrow(obs_keys),   "\n")
  cat("Panel not in Obs: ", nrow(p_not_o), "\n")
  cat("Obs not in Panel: ", nrow(o_not_p), "\n")
  if (nrow(p_not_o) > 0) {
    cat("\nExamples (Panel not in Obs):\n"); print(head(p_not_o, 8))
  }
  if (nrow(o_not_p) > 0) {
    cat("\nExamples (Obs not in Panel):\n"); print(head(o_not_p, 8))
  }

  # --- Merge (keep panel universe) ---
  merged <- panel %>%
    left_join(obs_raw, by = key) %>%
    select(-matches("_adm2_avg$"))  # just in case

  # Normalize any product tags that appear in names so they’re consistent
  # e.g., "era5-025" or "era5_025" -> "ERA5_025"
  nm <- names(merged)
  nm <- gsub("[-]", "_", nm)
  # Replace any variant of the product tag to canonical uppercase form
  # (safe no-op if not present)
  pat <- paste0("(?i)", gsub("_", "[-_]?", product_tag))
  nm <- gsub(pat, product_tag, nm, perl = TRUE)
  names(merged) <- nm

  # --- Post-merge diagnostics (generic; tavg_poly_1/lr_tavg names may vary) ---
  has_lr <- grep(paste0("(^|_)lr_tavg(_|$)"), names(merged), value = TRUE)
  has_t1 <- grep(paste0("(^|_)tavg_poly_1(_|$)"), names(merged), value = TRUE)

  diag <- merged %>%
    summarise(
      total_rows    = n(),
      with_lr_tavg  = if (length(has_lr) > 0) sum(!is.na(merged[[has_lr[1]]])) else NA_integer_,
      pct_lr_tavg   = if (length(has_lr) > 0) mean(!is.na(merged[[has_lr[1]]]))*100 else NA_real_,
      with_tpoly    = if (length(has_t1) > 0) sum(!is.na(merged[[has_t1[1]]])) else NA_integer_,
      pct_tpoly     = if (length(has_t1) > 0) mean(!is.na(merged[[has_t1[1]]]))*100 else NA_real_
    )
  cat("\n=== MERGE DIAG (", product_tag, ") ===\n", sep = ""); print(diag)

  # --- Write out ---
  out_path <- file.path(out_dir, paste0("global_mortality_panel_public_", product_tag, ".dta"))
  if (file.exists(out_path)) file.remove(out_path)
  write_dta(merged, out_path)
  message("✅ Wrote merged panel: ", out_path)
}

library(haven); library(dplyr)

df <- read_dta("/user/ab5405/summeraliaclimate/code/regressions/prep_panels/global_mortality_panel_public_MERRA2.dta")
t1 <- grep("(^|_)tavg_poly_1(_|$)", names(df), value=TRUE)[1]

miss <- df %>% 
  filter(is.na(.data[[t1]])) %>% 
  select(iso, adm1_id, adm2_id, year)

miss %>% count(iso, sort=TRUE) %>% print(n=99)       # where is it concentrated?
miss %>% count(iso, adm1_id, sort=TRUE) %>% head(30) # a few ADM1s often explain most

message("\nAll products processed.")
