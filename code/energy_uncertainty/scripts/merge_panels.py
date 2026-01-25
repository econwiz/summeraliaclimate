#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-

"""
merge_energy_with_climate.py

Attach product-specific country–year climate variables to the final
GMFD energy panel from the CIL release.

Inputs
------
1) Energy panel (GMFD-based), with country-year structure:
    /user/ab5405/summeraliaclimate/code/energy_consumption/\
        energy_data_release_2021oct21/DATA/regression/GMFD_TINV_clim_regsort.dta

2) Climate panels (built in Python) for each product:
    PROJECT_ROOT/data/country_climate/{PRODUCT}/
        {PRODUCT}_country_climate_1971_2010.csv

   where PRODUCT ∈ {GMFD, ERA5, JRA_3Q, MERRA2}

Outputs
-------
PROJECT_ROOT/data/merged_panels/{PRODUCT}_energy_panel.csv
PROJECT_ROOT/data/merged_panels/{PRODUCT}_energy_panel.dta
"""

from pathlib import Path
import pandas as pd

YEAR_MIN, YEAR_MAX = 1971, 2010

PROJECT_ROOT = Path("/user/ab5405/summeraliaclimate/code/energy_uncertainty")

# Where the new climate CSVs live
CLIMATE_DIR = PROJECT_ROOT / "data" / "country_climate"

# Where we write merged panels
MERGED_DIR = PROJECT_ROOT / "data" / "merged_panels"

# Final GMFD panel from the original CIL release
ENERGY_PANEL_PATH = Path(
    "/user/ab5405/summeraliaclimate/code/energy_consumption/"
    "energy_data_release_2021oct21/DATA/regression/GMFD_TINV_clim_regsort.dta"
)

PRODUCTS = ["GMFD", "ERA5", "JRA_3Q", "MERRA2"]


def main() -> None:
    print(f"[INFO] Reading base energy panel:\n  {ENERGY_PANEL_PATH}")
    energy = pd.read_stata(ENERGY_PANEL_PATH)

    if "country" not in energy.columns or "year" not in energy.columns:
        raise ValueError(
            "Expected 'country' and 'year' in energy panel. "
            f"Got: {list(energy.columns)}"
        )

    # Standardize merge keys
    energy = energy.copy()
    energy["iso"] = energy["country"]
    energy["year"] = energy["year"].astype(int)

    print("[INFO] Energy panel shape:", energy.shape)
    print(energy[["iso", "year"]].head())

    MERGED_DIR.mkdir(parents=True, exist_ok=True)

    for prod in PRODUCTS:
        print("\n" + "=" * 70)
        print(f"[INFO] Merging energy panel with {prod} climate")

        clim_dir = CLIMATE_DIR / prod
        clim_path = clim_dir / f"{prod}_country_year_{YEAR_MIN}_{YEAR_MAX}.csv"

        if not clim_path.exists():
            raise FileNotFoundError(f"Missing climate file for {prod}: {clim_path}")

        print(f"[INFO] Reading climate panel:\n  {clim_path}")
        climate = pd.read_csv(clim_path)

        # Basic checks
        if not {"iso", "year"}.issubset(climate.columns):
            raise ValueError(
                f"{prod} climate file must contain 'iso' and 'year'. "
                f"Columns are: {list(climate.columns)}"
            )

        print(f"[INFO] {prod} climate shape: {climate.shape}")

        # All climate columns to be attached (including TINV)
        climate_cols = [c for c in climate.columns if c not in ("iso", "year")]

        # Small debug: show TINV columns (if any)
        tinv_cols = [c for c in climate_cols if "TINV" in c]
        if tinv_cols:
            print(f"[INFO] {prod} TINV columns found:", tinv_cols)

        # Drop any existing columns with the same names to avoid duplicates
        energy_stripped = energy.drop(columns=climate_cols, errors="ignore")

        # Merge on iso-year
        merged = energy_stripped.merge(climate, on=["iso", "year"], how="left")

        n_total = len(merged)
        prod_cols = [c for c in merged.columns if c.endswith(f"_{prod}")]
        if prod_cols:
            n_missing = merged[prod_cols].isna().all(axis=1).sum()
            print(f"[INFO] {prod}: merged shape = {merged.shape}")
            print(
                f"[INFO] {prod}: rows with ALL {prod}-specific climate cols NA: "
                f"{n_missing} / {n_total}"
            )
        else:
            print(f"[WARN] No columns ending in _{prod} found after merge.")

        out_csv = MERGED_DIR / f"{prod}_energy_panel.csv"
        out_dta = MERGED_DIR / f"{prod}_energy_panel.dta"

        merged.to_csv(out_csv, index=False)
        print(f"[DONE] Wrote {out_csv}")

        merged.to_stata(out_dta, write_index=False)
        print(f"[DONE] Wrote {out_dta}")

    print("\n[INFO] Completed all merges.")


if __name__ == "__main__":
    main()
