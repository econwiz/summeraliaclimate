# Figure: Energy–Temperature Response by Income and Climate Tercile

Replicates Figure 1C from Rode et al. (2021), extended to four
observational climate products.

## What the figure shows

A 3×3 matrix of predicted energy–temperature response curves for
electricity demand. Rows separate income terciles (increasing income
from bottom to top); columns separate climate terciles by annual
cooling degree days (increasingly warm climate from left to right).
Each cell traces the predicted change in per-capita electricity
consumption (ΔE) relative to 20°C across the temperature range
−5°C to 35°C, evaluated at Year=2099 to capture the full temporal
trend.

## Script

`plot_energy_response.py`

## Inputs

| File | Source |
|------|--------|
| `data/regression/sters/FD_FGLS_inter_TINV_clim_quadinter_{PRODUCT}_coeff.csv` | Step 4 of pipeline |
| `output_original/sters/FD_FGLS_inter_TINV_clim_quadinter_coeff.csv` | CIL Zenodo release |
| `energy_data_release_2021oct21/DATA/regression/break_data_TINV_clim.dta` | CIL Zenodo release |

## Outputs (saved to `figures/comparison/`)

| File | Contents |
|------|----------|
| `NEW_all_products_electricity_TINV_clim_quadinter.pdf` | All four products overlaid |
| `ORIGINAL_electricity_TINV_clim_quadinter.pdf` | Original CIL estimates |
| `COMPARISON_GMFD_electricity_TINV_clim_quadinter.pdf` | Original vs new GMFD |

## Notes

- Income knot location (`ibar`) and tercile covariate means are read
  from `break_data_TINV_clim.dta`, matching the original CIL plotting
  script exactly.
- Original CIL coefficients use `_GMFD`-suffixed variable names; new
  pipeline coefficients are unsuffixed. The script detects this
  automatically.
- The `quadinter` spec includes temp×year and temp×year² interactions.
  If those terms were dropped as collinear during estimation, they
  contribute 0 to the plotted curve — this is correct behavior.