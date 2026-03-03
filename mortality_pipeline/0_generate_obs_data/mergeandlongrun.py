#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-
"""
Script B — Merge per-year outputs + ADM1 Long-Run Mean (Kevin-aligned)

- Reads yearly parquet outputs from Script A
- Concatenates into one DataFrame (ADM2-by-year)
- Computes ADM1 long-run mean temperature (1981–2010, or all available years)
- Fixes invalid geometries before dissolving
- Merges ADM1 long-run mean into panel universe
- Writes final wide CSV
"""

from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import xagg as xa

# ─── CONFIG ────────────────────────────────────────────────────────────────
PRODUCT       = "MERRA2"   # change as needed
TAG           = PRODUCT.replace("-", "_").upper()

BASE          = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data")
CAR_PATHS_CSV = BASE / "car_paths.csv"

PANEL_DTA     = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
SHAPEFILE     = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

IN_YRLY       = BASE / "obs_yearly"
OUTDIR        = BASE / "obs_csvs"; OUTDIR.mkdir(exist_ok=True)
OUT_CSV       = OUTDIR / f"{PRODUCT.replace('-', '_')}_by_region_year.csv"

# EU collapse (note Kevin excludes FRA but includes TUR)
EU_ISO = {
    "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST","FIN","GBR","GRC",
    "HRV","HUN","IRL","ISL","ITA","LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD",
    "NOR","POL","PRT","ROU","SVK","SVN","SWE","TUR"
}
def collapse_eu(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().where(~s.isin(EU_ISO), "EU")

# ─── 1. Mortality panel keys ──────────────────────────────────────────────
df_dta = pd.read_stata(PANEL_DTA)
years_insample = np.unique(df_dta.year.values)

# EU collapse
df_dta = df_dta.rename(columns={"iso":"iso0"})
df_dta["iso"] = df_dta.iso0.where(~df_dta.iso0.isin(EU_ISO), "EU")

keys  = df_dta[["iso","adm0","adm1_id","adm2_id","year"]].drop_duplicates()
years = np.sort(keys["year"].unique())

# ─── 2. Shapefile (insample regions) ──────────────────────────────────────
gdf = gpd.read_file(SHAPEFILE)
gdf["geometry"] = gdf["geometry"].make_valid()
gdf = gdf.to_crs("EPSG:4326")

# Subset to panel universe
gdf_sub = pd.merge(
    df_dta.loc[:,["iso","adm0","adm1_id","adm2_id"]].drop_duplicates(),
    gdf,
    on=["iso","adm1_id","adm2_id"]
)
gdf_sub = gdf_sub.set_geometry("geometry")

# Build dissolved ADM1
gdf_adm1 = gdf_sub.dissolve(by=["iso","adm1_id"]).reset_index().drop(columns="adm2_id")

# ─── 3. Load yearly parquet outputs ───────────────────────────────────────
files = sorted(IN_YRLY.glob(f"{PRODUCT.replace('-','_')}_*.parquet"))
if not files:
    raise FileNotFoundError(f"No yearly parquet files found in {IN_YRLY}")
dfs = [pd.read_parquet(p) for p in files]
df_by_year = pd.concat(dfs, ignore_index=True)
print(f"[INFO] Loaded {len(files)} yearly files → {len(df_by_year):,} rows")

# ─── 3. ADM1 Long-run Mean ─────────────────────────────────────────────────
paths = pd.read_csv(CAR_PATHS_CSV, dtype=str)
row   = paths.loc[paths["product"].str.strip() == PRODUCT].iloc[0]
tas_fp = row["tas_filepath"]

# Open tas and convert to °C
ds = xr.open_dataset(tas_fp, chunks={"lat": -1, "lon": -1, "time": 20})
ds["tas"] = ds["tas"].astype("float32") - 273.15

# Subset to long-run years (1981–2010, or all available)
years_lr = [y for y in years if 1981 <= y <= 2010] or list(years)
ds_lr_r = ds.sel(time=ds.time.dt.year.isin(years_lr))

# Compute annual means, then long-run mean
tmean_y  = ds_lr_r["tas"].groupby("time.year").mean("time")
tmean_lr = tmean_y.mean("year").rename("Tmean")

# Wrap into dataset so xagg can parse lat/lon correctly
tmean_lr = tmean_lr.to_dataset()

# Read shapefile and prep ADM1 polygons
gdf = gpd.read_file(SHAPEFILE).to_crs("EPSG:4326")
gdf["iso"] = collapse_eu(gdf["iso"])
gdf = gdf[gdf.geometry.notnull()]
gdf_invalid = gdf[~gdf.is_valid]
if not gdf_invalid.empty:
    print(f"[WARN] Found {len(gdf_invalid)} invalid geometries — fixing...")
    gdf.loc[gdf_invalid.index, "geometry"] = gdf_invalid.geometry.make_valid()

gdf_adm1 = gdf.dissolve(by=["iso", "adm1_id"], as_index=False)

# Compute weightmap + aggregation with xagg
with xa.set_options(silent=True):
    WM_ADM1 = xa.pixel_overlaps(tmean_lr, gdf_adm1)
with xa.set_options(impl="numba", silent=True):
    agg_lr = xa.aggregate(tmean_lr, WM_ADM1)

# Convert to DataFrame
df_lr = agg_lr.to_dataframe().reset_index()[["iso", "adm1_id", "Tmean"]]
df_lr = df_lr.rename(columns={"Tmean": f"lr_tavg_{TAG}_adm1_avg"})

# ─── 5. Merge and write ──────────────────────────────────────────────────
df_by_year["iso"]    = collapse_eu(df_by_year["iso"])
df_by_year["adm1_id"]= df_by_year["adm1_id"].astype(str)
df_by_year["adm2_id"]= df_by_year["adm2_id"].astype(str)

out = df_by_year.merge(df_lr, on=["iso","adm1_id"], how="left")
universe = keys[["iso","adm1_id","adm2_id","year"]].drop_duplicates()
out = universe.merge(out, on=["iso","adm1_id","adm2_id","year"], how="left")

out.to_csv(OUT_CSV, index=False)
print(f"✅ Wrote {OUT_CSV.name} with {len(out):,} rows → {OUTDIR}")
