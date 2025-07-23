# regressions/

This directory contains all of the Stata workflows for estimating
temperature–mortality response functions.

- **prep_data/**  
  Prepares the global mortality panel with ERA5, JRA‑3Q, MERRA2, GMFD.

- **age_combined_regressions/**  
  Runs 5 model specifications (pooled, FGLS, 13‑month) and saves `.ster` files.

- **output/**  
  Holds `.ster` regression output (ignored via `.gitignore`).

- **summarize_data.do**  
  Aggregates and compares regression summaries into tables/figures.

