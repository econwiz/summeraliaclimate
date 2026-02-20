# Pipeline Overview

## Step 1: Build country–year climate panel
**Script:** `/scripts/build_all_prod_regs.py`

For each product (ERA5, GMFD, JRA-3Q, MERRA2):
- Uses daily near-surface air temperature (tas) to compute annual
  temperature polynomials (temp1–4), HDD20, CDD20, polyAbove1–2,
  polyBelow1–2, and TINV cross terms (polyAbove×CDD, polyBelow×HDD).
  TINV cross terms are computed at the pixel level before aggregation —
  this is critical and cannot be done at the country level.
- Uses daily precipitation (pr, kg m⁻² s⁻¹), converted to mm/day
  (×86400), to compute annual total precipitation (precip1) and annual
  total squared daily precipitation (precip2).
- Aggregates all variables to country–year using GPW v4 year-2000
  population weights (static weights to avoid demographic confounding).

**Output:** `data/country_climate/{PRODUCT}/{PRODUCT}_country_year_1971_2010.csv`

## Step 2: Merge climate with energy panel
**Script:** `scripts/merge_panels.py`

Merges each product's country–year climate CSV into the CIL energy
panel on (iso, year) via left join, preserving the full energy panel
universe. For GMFD, freshly recomputed climate variables replace the
original CIL versions.

**Output:** one merged panel per product, ready for FD construction.

## Step 3: Construct first differences
**Script:** `build_FD_dataset_all.do` (Stata)

Generalizes the original CIL FD construction to all four products.
Computes first differences of load_pc, temperature polynomials,
precipitation, income, and all interaction terms (temp×year,
temp×year², temp×income spline, TINV cross terms). Outputs a
regression-ready panel per product.

**Output:** `data/regression/{PRODUCT}_TINV_clim_regsort.dta`

## Step 4: Run regressions and export results
**Scripts:** `FD_uninteracted_all.do`, `FD_interacted_all.do` (Stata)

For each product, estimates two specifications:
- **Uninteracted** (`2_FD_uninteracted_all.do`): quartic polynomial in
  temperature, no climate/income/year interactions. Two-stage FGLS
  with population-weighted first stage.
- **Interacted** (`3_FD_interacted_all.do`): full CIL specification
  with TINV climate heterogeneity, income spline, and year² interactions
  (quadinter). Two-stage FGLS with unweighted first stage, matching
  the original CIL estimation procedure.

Each script saves Stata `.ster` files and exports coefficient vectors
and variance-covariance matrices to `.csv` for use in Python plotting.

**Output:** `data/regression/sters/FD_FGLS_inter_{MODEL}_{PRODUCT}_coeff.csv`
            `data/regression/sters/FD_FGLS_inter_{MODEL}_{PRODUCT}_vcov_long.csv`