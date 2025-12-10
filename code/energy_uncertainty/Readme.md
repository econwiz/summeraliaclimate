# Overview

1. **Build country–year climate**
   - Script: `build_PRODUCT_country_climate_multi.py`
   - For each product (ERA5-025, GMFD, JRA-3Q, MERRA2):
     - Use **daily tas** to build annual temp polynomials, HDD20, CDD20, polyAbove*, polyBelow*.
     - Use **monthly pr_Amon (kg m-2 s-1)**, convert to **mm/month**, take annual mean ?? this is the confusing part.
     - Aggregate all variables to **country–year** with GPW population weights.
     - Output: `data/country_climate/{PRODUCT}/{PRODUCT}_country_year_1971_2010.csv`.

2. **Merge with energy panel**
   - Script: climate–energy merge script (Python).
   - Merge each `{PRODUCT}_country_year_*.csv` into the GMFD energy panel on `(iso, year)`.
   - Output: combined panel with energy + climate variables for all products.

3. **Construct first differences (FD)**
   - Script: `build_FD_dataset_all.do` (Stata).
   - This is from the original pipeline to get a regression-ready panel.

4. **Run regressions**
   - Scripts:  `uninteractedreg.do` and `interactedreg.do`
   - Save results as `.ster` and export coefficient tables to `.csv`?? not sure if this is correct either.
