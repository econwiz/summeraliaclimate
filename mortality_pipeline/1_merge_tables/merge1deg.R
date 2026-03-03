#!/usr/bin/env Rscript

# FLEXIBLE RESOLUTION VERSION
# 
# Merges ADM2-level climate observations (from 0_generate_obs_data) into the
# Carleton et al. global mortality panel. Supports both 0.25° and 1° resolutions.
#
# CHANGE LOG:
# - Added RESOLUTION parameter to control input/output directories
# - Outputs go to prep_panels/1deg/ or prep_panels/025deg/
# - Auto-detects all paroducts in the resolution-specific directory
#
# For each *_by_region_year.csv found in obs_dir, the script:
#   1. Normalizes join keys (ISO, ADM1/ADM2 IDs, year) in both the panel and obs data
#   2. Runs pre-merge diagnostics to flag unmatched keys
#   3. Left-joins obs onto the panel (preserving the full panel universe)
#   4. Standardizes product tag formatting in column names
#   5. Runs post-merge diagnostics to check coverage
#   6. Writes the merged panel as a .dta file to out_dir
#
# Inputs:  global_mortality_panel_public.dta, <product>_by_region_year.csv
# Outputs: global_mortality_panel_public_<PRODUCT>.dta

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(haven)
  library(janitor); library(stringr); library(tools)
})

# ========================================================================
# USER SETTINGS
# ========================================================================

panel_in <- "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"

# RESOLUTION: "025deg" or "1deg"
RESOLUTION <- "1deg"  # ← CHANGE THIS to switch between resolutions

# Input directory (resolution-specific)
obs_dir  <- paste0("/user/ab5405/summeraliaclimate/code/mortality_pipeline/0_generate_obs_data/obs_csvs/", RESOLUTION)

# Output directory (resolution-specific)
out_dir  <- paste0("/user/ab5405/summeraliaclimate/code/mortality_pipeline/prep_panels/", RESOLUTION)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", strrep("=", 70), "\n")
cat("CONFIGURATION:\n")
cat("  Resolution:", RESOLUTION, "\n")
cat("  Input dir:", obs_dir, "\n")
cat("  Output dir:", out_dir, "\n")
cat(strrep("=", 70), "\n\n")

# If you want to restrict to specific products, put them here
# (otherwise auto-detects all products in the directory):
only_products <- c()  # Empty = process all products
# Examples: c("ERA5_025","GMFD","JRA_3Q","MERRA2")

# ========================================================================
# EU COLLAPSE
# ========================================================================

euro_iso <- c(
  "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP",
  "EST","FRA","FIN","GBR","GRC","HRV","HUN","IRL","ISL","ITA",
  "LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
  "POL","PRT","ROU","SVK","SVN","SWE","TUR"
)

# ========================================================================
# LOAD PANEL
# ========================================================================

panel <- read_dta(panel_in) %>%
  clean_names() %>%
  mutate(
    iso     = toupper(str_trim(as.character(iso))),
    iso     = if_else(iso %in% euro_iso, "EU", iso),
    adm1_id = str_trim(as.character(adm1_id)),
    adm2_id = str_trim(as.character(adm2_id)),
    year    = as.integer(year)
  ) %>%
  # Drop any pre-existing climate columns
  select(-matches("era5[-_]?025|merra2|gmfd|jra|best", ignore.case = TRUE))

# ========================================================================
# FIND OBS CSVs
# ========================================================================

obs_files <- list.files(obs_dir, pattern = "_by_region_year\\.csv$", full.names = TRUE)

if (length(only_products) > 0) {
  keep <- vapply(obs_files, function(p) {
    base <- basename(p)
    prod <- toupper(gsub("-", "_", sub("_by_region_year\\.csv$", "", base)))
    prod %in% toupper(gsub("-", "_", only_products))
  }, logical(1))
  obs_files <- obs_files[keep]
}

if (length(obs_files) == 0) {
  stop("No observation files found in: ", obs_dir)
}

cat("Found", length(obs_files), "product(s) to merge:\n")
for (f in obs_files) cat("  -", basename(f), "\n")
cat("\n")

# ========================================================================
# HELPER FUNCTION
# ========================================================================

infer_product <- function(path) {
  nm <- basename(path)
  toupper(gsub("-", "_", sub("_by_region_year\\.csv$", "", nm)))
}

# ========================================================================
# MERGE LOOP
# ========================================================================

key <- c("iso","adm1_id","adm2_id","year")

for (obs_csv in obs_files) {
  product_tag <- infer_product(obs_csv)
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

  # --- Pre-merge diagnostics ---
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

  # --- Merge ---
  merged <- panel %>%
    left_join(obs_raw, by = key) %>%
    select(-matches("_adm2_avg$"))

  # Normalize product tags in column names
  nm <- names(merged)
  nm <- gsub("[-]", "_", nm)
  pat <- paste0("(?i)", gsub("_", "[-_]?", product_tag))
  nm <- gsub(pat, product_tag, nm, perl = TRUE)
  names(merged) <- nm

  # --- Post-merge diagnostics ---
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

  # --- Write output ---
  out_path <- file.path(out_dir, paste0("global_mortality_panel_public_", product_tag, ".dta"))
  if (file.exists(out_path)) file.remove(out_path)
  write_dta(merged, out_path)
  message("✅ Wrote merged panel: ", out_path)
}

cat("\n", strrep("=", 70), "\n")
cat("All products processed successfully!\n")
cat("Resolution:", RESOLUTION, "\n")
cat("Output location:", out_dir, "\n")
cat(strrep("=", 70), "\n\n")