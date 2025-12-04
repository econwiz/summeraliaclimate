#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-

"""
merge_energy_with_climate.py

Goal
----
Take the FINAL GMFD energy panel from the CIL repo and attach
climate regressors from each product-specific country-year panel:

    GMFD_country_climate_1971_2010.csv
    ERA5_country_climate_1971_2010.csv
    JRA_3Q_country_climate_1971_2010.csv
    (MERRA2_country_climate_1971_2010.csv — later)

Output
------
energy_uncertainty/data/merged_panels/{PRODUCT}_energy_panel.csv

where PRODUCT ∈ {GMFD, ERA5, JRA_3Q, (MERRA2)}
"""

from pathlib import Path
import pandas as pd
import os

YEAR_MIN, YEAR_MAX = 1971, 2010

PROJECT_ROOT = Path("/user/ab5405/summeraliaclimate/code/energy_uncertainty")

CLIMATE_DIR = PROJECT_ROOT / "data" / "country_year"

MERGED_DIR = PROJECT_ROOT / "data" / "merged_panels"

ENERGY_PANEL_PATH = Path("/user/ab5405/summeraliaclimate/code/energy_consumption/energy_data_release_2021oct21/DATA/regression/GMFD_TINV_clim_regsort.dta")

PRODUCTS = [
    "GMFD",
    "ERA5",
    "JRA_3Q",
    # "MERRA2",  # ← uncomment once MERRA2 climate CSV exists
]

def main():
    print(f"[INFO] Reading energy panel from:\n {ENERGY_PANEL_PATH}")
    energy = pd.read_stata(ENERGY_PANEL_PATH)

    if "country" not in energy.columns or "year" not in energy.columns:
        raise ValueError(
            f"Expected 'country' and 'year' in energy panel. Got: {energy.columns.tolist()}"
        )

    energy["iso"] = energy["country"]
    energy["year"] = energy["year"].astype(int)

    print("[INFO] Energy panel shape:", energy.shape)
    print("[INFO] Example rows:")
    print(energy[["iso", "year"]].head())

    for prod in PRODUCTS:
        print("\n" + "=" * 70)
        print(f"[INFO] Merging energy with {prod} climate panel")

        clim_dir = CLIMATE_DIR / prod
        clim_path = clim_dir / f"{prod}_country_climate_{YEAR_MIN}_{YEAR_MAX}.csv"

        print(f"[INFO] Reading climate panel from:\n {clim_path}")
        climate = pd.read_csv(clim_path)

        print("[INFO]", prod, "climate shape:", climate.shape)

        climate_cols = [c for c in climate.columns if c not in ("iso", "year")]

        # Drop any of those climate columns from energy if they already exist
        energy_stripped = energy.drop(columns=climate_cols, errors="ignore")
        
        merged = energy_stripped.merge(climate, on=["iso", "year"], how = "left")

        n_total = len(merged)
        n_missing = merged[merged.filter(like=prod).isna().all(axis=1)].shape[0]

        print(f"[INFO] {prod}: merged shape = {merged.shape}")
        print(f"[INFO] {prod}: rows with ALL {prod} climate cols NA: {n_missing} / {n_total}")

        out_csv = MERGED_DIR / f"{prod}_energy_panel.csv"
        out_dta = MERGED_DIR / f"{prod}_energy_panel.dta"
        
        merged.to_csv(out_csv, index=False)
        print(f"[DONE] Wrote {out_csv}")
        
        # ALSO write Stata .dta to avoid CSV/BOM header issues
        merged.to_stata(out_dta, write_index=False)
        print(f"[DONE] Wrote {out_dta}")


    print("\n[All Done]")

if __name__ == "__main__":
    main()
        