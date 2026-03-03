#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8

import xarray as xr
import pandas as pd
from pathlib import Path
import numpy as np
import os

BASE = Path(__file__).resolve().parent
paths = pd.read_csv(BASE / "car_paths.csv")

def is_missing(x):
    return x is None or (isinstance(x, float) and np.isnan(x)) or str(x).strip() == ""

def to_path_str(x):
    return str(x).strip()

def open_climate(path_str):
    if path_str.endswith(".zarr"):
        return xr.open_zarr(path_str, consolidated=False)
    return xr.open_dataset(path_str)

for idx, row in paths.iterrows():
    prod = str(row.get("product", "")).strip() or f"(row {idx})"

    tas_raw = row.get("tas_filepath", None)
    pr_raw  = row.get("precip_filepath", None)

    # Always check tas
    if is_missing(tas_raw):
        print(f"\n[{prod}]")
        print("  SKIP: missing tas_filepath")
        continue

    tas_path = to_path_str(tas_raw)

    if not (tas_path.endswith(".zarr") or os.path.exists(tas_path)):
        print(f"\n[{prod}]")
        print("  SKIP: tas file not found:", tas_path)
        continue

    try:
        ds_t = open_climate(tas_path)
    except Exception as e:
        print(f"\n[{prod}]")
        print("  ERROR opening tas dataset:", e)
        continue

    print(f"\n[{prod}]")
    try:
        print("  tas dtype:", ds_t["tas"].dtype)
    except Exception as e:
        print("  tas dtype: ERROR:", e)

    # Skip precip for MERRA2, otherwise check
    if prod.upper().startswith("MERRA"):
        print("  precip: SKIPPED (MERRA has no precip data)")
        continue

    if is_missing(pr_raw):
        print("  SKIP: missing precip_filepath")
        continue

    pr_path = to_path_str(pr_raw)
    if not (pr_path.endswith(".zarr") or os.path.exists(pr_path)):
        print("  SKIP: precip file not found:", pr_path)
        continue

    try:
        ds_p = open_climate(pr_path)
    except Exception as e:
        print("  ERROR opening precip dataset:", e)
        continue

    pvar = next((v for v in ("prcp","pr","precip","tp") if v in ds_p.data_vars), None)
    if pvar:
        print(f"  {pvar} dtype:", ds_p[pvar].dtype)
    else:
        print("  precip var: NOT FOUND (vars: {})".format(list(ds_p.data_vars)))
