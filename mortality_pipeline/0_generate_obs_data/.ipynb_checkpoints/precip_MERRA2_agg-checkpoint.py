#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8

import xagg as xa
import geopandas as gpd
import xarray as xr
import pandas as pd
from pathlib import Path

product = "MERRA2"
BASE = Path(__file__).resolve().parent

# Load paths
df_paths = pd.read_csv(BASE / "car_paths.csv", dtype=str)
df_paths["product"] = df_paths["product"].str.strip()
row_merra = df_paths.loc[df_paths["product"] == "MERRA2"].iloc[0]
row_gmfd  = df_paths.loc[df_paths["product"] == "GMFD"].iloc[0]

# Load shapefile and fix geometries
shape_path = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
gdf_adm2 = gpd.read_file(shape_path).to_crs("EPSG:4326")
gdf_adm2["geometry"] = gdf_adm2.buffer(0)
gdf_adm2 = gdf_adm2[gdf_adm2.is_valid & ~gdf_adm2.is_empty]
gdf_adm1 = gdf_adm2.dissolve(by=["iso", "adm1_id"], as_index=False)

def open_climate(path):
    return xr.open_zarr(path, consolidated=False) if path.endswith(".zarr") else xr.open_dataset(path)

# Load climate data with optimized chunking
ds_tas = open_climate(row_merra["tas_filepath"]).chunk({"time": 365})
ds_tas = ds_tas.assign(tas=ds_tas.tas - 273.15)

ds_pr  = open_climate(row_gmfd["precip_filepath"]).chunk({"time": 365})
precip_var = next(v for v in ["prcp","pr","precip"] if v in ds_pr.data_vars)

# Precompute weightmaps
sample_tas = ds_tas.isel(time=0, drop=True)
sample_pr  = ds_pr.isel(time=0, drop=True)

wm_adm2_tas = xa.pixel_overlaps(sample_tas, gdf_adm2)
wm_adm1_tas = xa.pixel_overlaps(sample_tas, gdf_adm1)
wm_adm2_pr  = xa.pixel_overlaps(sample_pr, gdf_adm2)

# ADM1 long-run average temp
tmean_lr = ds_tas.sel(time=slice("1981-01-01", "2010-12-31")).tas.mean("time")
with xa.set_options(impl="numba", silent=True):
    agg_lr_tas = xa.aggregate(tmean_lr, wm_adm1_tas)

df_lr_tas = (agg_lr_tas
             .to_dataframe().reset_index()
             .rename(columns={'tas': "lr_tavg_MERRA2_adm1_avg"}))

years = range(1981, 2014)
OUTDIR = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs")
out_csv = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"

# Aggregate MERRA2 Temperature
rows_tas = []
for yr in years:
    ds_yr_tas = ds_tas.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))
    ds_poly_tas = xr.Dataset({
        "T1_sum": ds_yr_tas.tas.sum("time"),
        "T2_sum": (ds_yr_tas.tas ** 2).sum("time"),
        "T3_sum": (ds_yr_tas.tas ** 3).sum("time"),
        "T4_sum": (ds_yr_tas.tas ** 4).sum("time"),
    })

    with xa.set_options(impl="numba", silent=True):
        agg_tas = xa.aggregate(ds_poly_tas, wm_adm2_tas)

    df_tas = agg_tas.to_dataframe().reset_index()
    df_tas['year'] = yr
    rows_tas.append(df_tas)

df_tas_full = pd.concat(rows_tas, ignore_index=True)
df_tas_full = df_tas_full.merge(df_lr_tas, on=['iso', 'adm1_id'], how='left')

# Aggregate GMFD Precipitation
rows_pr = []
for yr in years:
    ds_yr_pr = ds_pr.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))
    pr_sum = ds_yr_pr[precip_var].sum("time")

    with xa.set_options(impl="numba", silent=True):
        agg_pr = xa.aggregate(pr_sum, wm_adm2_pr)

    df_pr = agg_pr.to_dataframe().reset_index().rename(columns={precip_var: "prcp_annual_GMFD"})
    df_pr['year'] = yr

    for p in [1, 2, 3, 4]:
        df_pr[f"prcp_poly_{p}_GMFD"] = df_pr["prcp_annual_GMFD"] ** p

    rows_pr.append(df_pr)

df_pr_full = pd.concat(rows_pr, ignore_index=True)

# Final merge of temp and precip data
df_final = pd.merge(df_tas_full, df_pr_full, on=['iso', 'adm1_id', 'adm2_id', 'year'], how='inner')

# Write final CSV
df_final.to_csv(out_csv, index=False)
print(f"✅ Wrote {out_csv.name} ({len(df_final)} rows)")
