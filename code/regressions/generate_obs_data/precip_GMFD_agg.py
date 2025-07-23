#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8

import numpy as np
import xagg as xa
import geopandas as gpd
import xarray as xr
import pandas as pd
from pathlib import Path

product = "GMFD"
BASE    = Path(__file__).resolve().parent

# your existing paths setup
df_paths   = pd.read_csv(BASE / "car_paths.csv", dtype=str)
df_paths["product"] = df_paths["product"].str.strip()
row        = df_paths.loc[df_paths["product"] == product].iloc[0]

shape_path = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
gdf        = gpd.read_file(shape_path, engine="fiona").to_crs("EPSG:4326")
gdf["geometry"] = gdf.buffer(0)
gdf = gdf[gdf.is_valid & ~gdf.is_empty]

def open_climate(path):
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path)

# 1) load temperature (and possibly precipitation) dataset
ds_full = open_climate(row["filepath"]).chunk({"time": 1})

# 2) inspect to find precip var name  <— NEW
print("Data variables available:", list(ds_full.data_vars))
# common names might be 'prcp' or 'pr' or 'precip'; adjust below if needed
precip_var = "prcp" if "prcp" in ds_full.data_vars else "pr"

# convert temperature to °C
ds_full = ds_full.assign(tas = ds_full.tas - 273.15)

# build weight map
sample = ds_full.isel(time=0, drop=True)
wm     = xa.pixel_overlaps(sample, gdf)

years = range(1981, 2014)
OUTDIR  = BASE / "obs_csvs_prceip"
OUTDIR.mkdir(parents=True, exist_ok=True)
out_csv = OUTDIR / f"{product.replace('-','_')}_by_region_year.csv"

rows = []
for yr in years:
    print(f"  Year {yr}")
    ds_yr = ds_full.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))

    # —— existing temperature metrics —— 
    t1    = ds_yr.tas.sum("time").rename("T1_sum")
    t2    = (ds_yr.tas**2).sum("time").rename("T2_sum")
    t3    = (ds_yr.tas**3).sum("time").rename("T3_sum")
    t4    = (ds_yr.tas**4).sum("time").rename("T4_sum")
    tmean = ds_yr.tas.mean("time").rename("Tmean")

    # —— NEW: annual-mean precipitation —— 
    pr_mean = ds_yr[precip_var].mean("time").rename("PRCP_mean")

    # merge into one xarray Dataset
    ds_poly = xr.merge([t1, t2, t3, t4, tmean, pr_mean]) \
                .chunk({"lat": -1, "lon": -1})

    # aggregate to ADM1
    with xa.set_options(impl="numba", silent=True):
        agg = xa.aggregate(ds_poly, wm)

    # to pandas and rename columns
    df_reg = (
        agg.to_dataframe()
           .reset_index()
           .rename(columns={
               "T1_sum":    f"tavg_poly_1_{product.upper()}",
               "T2_sum":    f"tavg_poly_2_{product.upper()}",
               "T3_sum":    f"tavg_poly_3_{product.upper()}",
               "T4_sum":    f"tavg_poly_4_{product.upper()}",
               "Tmean":     f"lr_tavg_{product.upper()}_adm1_avg",
               "PRCP_mean": f"prcp_annual_{product.upper()}"     # <— NEW
           })
           .assign(year=yr, product=product)
    )

    # —— NEW: build precip polynomials in pandas ——
    for p in [1,2,3,4]:
        df_reg[f"prcp_poly_{p}_{product.upper()}"] = df_reg[f"prcp_annual_{product.upper()}"]**p

    rows.append(df_reg)

# finalize
prod_df = pd.concat(rows, ignore_index=True).fillna(0)
print("Writing output to:", out_csv)
prod_df.to_csv(out_csv, index=False)
print(f"↳ Wrote {out_csv.name} ({len(prod_df)} rows)")
print("✅ Done!")
