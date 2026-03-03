# regressions/

This directory contains all of the Stata workflows for estimating
temperature–mortality response functions.

- **0_generate_obs_data/**  
  Generates observational climate variables aggregated to ADM2/ADM1.  
  Climate products and file paths are defined in `car_paths.csv` in this directory.

- **1_merge_tables/**  
  Merges generated climate data with mortality and population panels at the ADM2 level, then collapses to ADM1 where needed.

- **2_prep_data/**  
  Cleans merged tables: creates exposure bins, handles missing values, and formats variables for regression.  
  Outputs intermediate Stata panels (`.dta`) for the merged data.

- **3_age_combined_regressions/**  
  Runs five model specifications (pooled, fixed effects, FGLS, 13‑month lag, cluster-robust variants) using combined age exposures.  
  Saves `.ster` and `.csvv` files in `output/`.

- **4_age_interacted_regressions/**  
  Runs regressions with interaction terms between age groups and exposure bins.  
  Saves `.ster` and `.csvv` files in `output/`.

- **prep_data/**  
  Contains the Stata panels (`.dta`) created immediately after the merge step, ready for regression scripts.

- **output_panels/**  
  Stores final regression-ready panels (`.dta`) after all prep steps.

- **output/**  
  Holds regression output files: `.ster`, `.csvv`, and log files.  
