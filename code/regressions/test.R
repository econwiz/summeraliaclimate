#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(haven); library(stringr); library(tidyr)
})

# --- paths (adjust product_tag if comparing a different product) ---
product_tag <- "JRA_3Q"
mine_path   <- paste0("/user/ab5405/summeraliaclimate/code/regressions/prep_panels/global_mortality_panel_public_", product_tag, ".dta")
kev_path    <- file.path("/shared/share_hle/data/aux_data", paste0("reginputs_ann_", gsub("_","-",product_tag), "_regvars.csv"))

stopifnot(file.exists(mine_path), file.exists(kev_path))

# --- helpers ---
euro_iso <- c("AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST","FRA","FIN","GBR","GRC",
              "HRV","HUN","IRL","ISL","ITA","LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
              "POL","PRT","ROU","SVK","SVN","SWE","TUR")
collapse_eu <- function(s) { s <- toupper(trimws(as.character(s))); ifelse(s %in% euro_iso, "EU", s) }
norm_ids <- function(df) {
  df %>% mutate(
    iso = collapse_eu(iso),
    adm1_id = trimws(as.character(adm1_id)),
    adm2_id = trimws(as.character(adm2_id)),
    year = as.integer(year)
  )
}
days_in_year <- function(y) ifelse((y %% 400 == 0) | (y %% 4 == 0 & y %% 100 != 0), 366, 365)

# --- load ---
mine_raw <- read_dta(mine_path)
# keep TOTAL safely: in Stata it’s a labelled numeric where 0 == "total"
if ("agegroup" %in% names(mine_raw) && inherits(mine_raw$agegroup, "haven_labelled")) {
  mine_raw <- dplyr::filter(mine_raw, as.numeric(agegroup) == 0)
} else if ("agegroup" %in% names(mine_raw)) {
  mine_raw <- dplyr::filter(mine_raw, tolower(trimws(as.character(agegroup))) == "total")
}
mine <- mine_raw %>% norm_ids()

kev <- read_csv(kev_path, show_col_types=FALSE) %>% norm_ids()

# --- restrict to Kevin’s key universe ---
key <- c("iso","adm1_id","adm2_id","year")
kev_keys  <- kev %>% distinct(across(all_of(key)))
mine_keys <- mine %>% distinct(across(all_of(key)))
mine <- semi_join(mine, kev_keys, by = key)

# --- map mine var names (they carry product tag) to Kevin’s (generic) ---
v_tavg <- paste0("tavg_poly_", 1:4)
v_prcp <- paste0("prcp_poly_", 1:4)

mine_cols <- c(
  setNames(paste0(v_tavg, "_", product_tag), v_tavg),
  setNames(paste0(v_prcp, "_", product_tag), v_prcp)
)

# lr_tavg column in “mine” often looks like: lr_tavg_<TAG>_adm1_avg (or similar)
mine_lr <- grep(paste0("^lr_tavg.*", product_tag, ".*adm1.*avg$|^lr_tavg.*", product_tag, "$"), names(mine), value=TRUE, ignore.case=TRUE)
if (length(mine_lr) == 0) {
  stop("Could not find lr_tavg column in mine for product ", product_tag)
}
mine_lr <- mine_lr[1]  # take the first match
kev_lr  <- "lr_tavg"

# --- join value tables on the key ---
# Before the inner_join in your checker:
mine_vals <- mine %>%
  distinct(iso, adm1_id, adm2_id, year, .keep_all = TRUE)

kev_vals  <- kev %>%
  distinct(iso, adm1_id, adm2_id, year, .keep_all = TRUE)

cmp <- inner_join(mine_vals, kev_vals, by = c("iso","adm1_id","adm2_id","year"),
                  suffix = c(".mine",".kev"))

# --- function to compare a pair of columns and print diagnostics ---
compare_var <- function(base_name, mine_name, kev_name, df) {
  if (!(mine_name %in% names(df)) || !(kev_name %in% names(df))) {
    cat("\n[", base_name, "] Skipping (missing one or both columns)\n", sep = "")
    return(invisible(NULL))
  }
  x <- df[[mine_name]]
  y <- df[[kev_name]]
  na_mismatch <- sum(is.na(x) != is.na(y))
  both <- complete.cases(x, y)
  n <- sum(both)
  if (n == 0) {
    cat("\n[", base_name, "] No overlapping non-NA rows; NA mismatch:", na_mismatch, "\n", sep = "")
    return(invisible(NULL))
  }
  xd <- x[both]; yd <- y[both]
  diff <- xd - yd
  rel  <- ifelse(abs(yd) > 0, diff / pmax(1e-12, yd), NA_real_)
  ratio <- ifelse(abs(yd) > 0, xd / yd, NA_real_)

  cat("\n=== ", base_name, " ===\n", sep = "")
  cat("rows:", n, " | NA mismatch:", na_mismatch, "\n")
  cat("abs diff  median:", signif(median(abs(diff), na.rm=TRUE),6),
      " max:", signif(max(abs(diff), na.rm=TRUE),6), "\n")
  cat("rel diff  median:", signif(median(abs(rel), na.rm=TRUE),6),
      " p99:", signif(quantile(abs(rel), 0.99, na.rm=TRUE),6), "\n")
  cat("corr(x_mine, y_kev):", signif(cor(xd, yd), 8), "\n")
  cat("ratio median (mine/kev):", signif(median(ratio, na.rm=TRUE),8), "\n")

  # special detectors
  if (grepl("^tavg_poly_1$", base_name)) {
    # if Kelvin vs Celsius: (x - y)/days ≈ 273.15
    yrs <- df$year[both]
    days <- days_in_year(yrs)
    off <- (diff / days)
    m_off <- median(off, na.rm=TRUE)
    if (is.finite(m_off) && abs(m_off - 273.15) < 0.5) {
      cat("⚠︎ Looks like K vs °C: median (mine - kev)/days ≈ 273.15\n")
    }
  }
  if (grepl("^prcp_poly_", base_name)) {
    med_ratio <- median(ratio, na.rm=TRUE)
    if (is.finite(med_ratio) && (abs(med_ratio - 86400) < 1 || abs(med_ratio - 1/86400) < 1e-4)) {
      cat("⚠︎ Looks like precip unit mismatch (kg m^-2 s^-1 vs mm/day): median mine/kev ≈ ", med_ratio, "\n", sep="")
    }
  }

  # print top 10 worst rows
  worst <- order(-abs(diff))[1:min(10, length(diff))]
  ex <- df[both, c("iso","adm1_id","adm2_id","year")]
  ex$mine <- xd; ex$kev <- yd; ex$diff <- diff
  cat("Worst 10 by |diff|:\n")
  print(ex[worst, ], n=10)
}

# --- run comparisons ---
for (k in 1:4) {
  compare_var(paste0("tavg_poly_", k),
              paste0("tavg_poly_", k, "_", product_tag),
              paste0("tavg_poly_", k),
              cmp)
}
for (k in 1:4) {
  compare_var(paste0("prcp_poly_", k),
              paste0("prcp_poly_", k, "_", product_tag),
              paste0("prcp_poly_", k),
              cmp)
}
compare_var("lr_tavg", mine_lr, kev_lr, cmp)
