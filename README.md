# summeraliaclimate

This repo bundles **two independent toolkits**:

| path | purpose | quick run |
|------|---------|-----------|
| `code/regressions/` | Re‑run (and extend) the temperature‑mortality regressions of **Carleton et al. 2022** with new observational data sets (ERA5‑0.25, GMFD, JRA‑3Q, MERRA2). | ```bash\nconda activate hle_iv\n# 1 build covariates for each product\npython code/regressions/0_generate_obs_data/precip_ERA5_agg.py  # repeat for other products\n# 2 merge obs → mortality panel\npython code/regressions/1_merge_tables/merge_era_new.R\n# 3 run age combined Stata regression specs 1‑5\nstata -b do code/regressions/3_age_combined_regressions/age_combined_regressions_era5.do\n# 4 run age spec interacted Stata regression spec 2\nstata -b do code/regressions/4_age_interacted_regressions/age_spec_interacted_regressions_era5.do\n``` |
| `code/var_partitioning/` | Split ensemble variance in any impact field into **model** vs **internal** components and plot bar / map figures. | ```bash\nconda activate hle_iv\ncd code/var_partitioning\npython generate_figures.py               # main driver\n# or individual scripts\python gen_abs_var_figs.py\npython gen_bar_graphs.py\n``` |

---

## Outputs

* **Regressions:**  
  * `.ster` binary outputs in `code/regressions/output/age_combined_*`  
  * '.csvv' outputs in 'code/regressions/output/age_combined_csvv' 'code/regressions/output/age_spec_interacted_csvv'/
* **Variance partitioning:**  
  * all PNG / PDF figures in `code/var_partitioning/figures/`

---

Questions?  
*Re‑run pipeline* → **Alia Bonanno** (`ab5405@columbia.edu`)  
*Variance partitioning* → open an issue or ping **@econwiz**

