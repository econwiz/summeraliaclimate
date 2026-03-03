#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-

"""
Climate Aggregation Script — ADM2-by-Year + ADM1 Long-Run Mean (Carleton-style)

- Handles GMFD, ERA5-025, etc.
- Aggregates tas and precip (power 1–4) over ADM2-by-year
- Computes ADM1 long-run mean temperature
- Panel-anchored output with EU collapse
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import xesmf as xe
import xagg as xa
from xagg.core import numba_aggregate
import dask

# ─── Dask setup ────────────────────────────────────────────────────────────────
dask.config.set({"array.slicing.split_large_chunks": True})
try:
    from pyproj import datadir
    os.environ["PROJ_LIB"] = datadir.get_data_dir()
except Exception:
    pass

# ─── User Settings ─────────────────────────────────────────────────────────────
PRODUCT       = "GMFD"
BASE = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data")
CAR_PATHS_CSV = BASE / "car_paths.csv"

PANEL_DTA     = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
SHAPEFILE     = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
OUTDIR        = BASE / "obs_csvs"
OUT_CSV       = OUTDIR / f"{PRODUCT.replace('-', '_')}_by_region_year.csv"
OUTDIR.mkdir(parents=True, exist_ok=True)

LAT_TARGET = np.arange(-56 + 0.125, 86 - 0.125, 0.25, dtype="float32")
LON_TARGET = np.arange(-180 + 0.125, 180 - 0.125, 0.25, dtype="float32")
OUT_CHUNKS = {"lat": 40, "lon": 40}

EU_ISO = {
    "AUT", "BEL", "BGR", "CHE", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FRA", "FIN",
    "GBR", "GRC", "HRV", "HUN", "IRL", "ISL", "ITA", "LIE", "LTU", "LUX", "LVA", "MKD",
    "MLT", "MNE", "NLD", "NOR", "POL", "PRT", "ROU", "SVK", "SVN", "SWE", "TUR"
}

def collapse_eu(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().where(~s.isin(EU_ISO), "EU")

def tag(name: str) -> str:
    return name.replace("-", "_").upper()

TAG = tag(PRODUCT)

def is_nonempty(x) -> bool:
    return bool(x) and str(x).strip().lower() not in {"", "none", "nan"}

def open_any(path: str) -> xr.Dataset:
    return xr.open_zarr(path, consolidated=False) if path.endswith(".zarr") else xr.open_dataset(path, chunks={"time": 30, "lat": 120, "lon": 120})

# ─── xESMF setup ───────────────────────────────────────────────────────────────
def _coerce_latlon(ds: xr.Dataset) -> xr.Dataset:
    if "latitude" in ds.coords: ds = ds.rename({"latitude": "lat"})
    if "longitude" in ds.coords: ds = ds.rename({"longitude": "lon"})
    if "lon" in ds:
        lon = ds["lon"].values
        if lon.min() >= 0 and lon.max() > 180:
            lon = ((lon + 180) % 360) - 180
            ds = ds.assign_coords(lon=("lon", lon))
        if np.any(np.diff(ds["lon"].values) < 0):
            ds = ds.sortby("lon")
    if "lat" in ds.coords and np.any(np.diff(ds["lat"].values) < 0):
        ds = ds.sortby("lat")
    return ds

def make_target_grid() -> xr.Dataset:
    return xr.Dataset(coords={"lat": LAT_TARGET, "lon": LON_TARGET})

def build_regridder(src_like: xr.Dataset) -> xe.Regridder:
    src = _coerce_latlon(src_like)
    return xe.Regridder(src, make_target_grid(), method="bilinear", periodic=True, ignore_degenerate=False)

def regrid_xesmf(ds_in: xr.Dataset, rgrd: xe.Regridder) -> xr.Dataset:
    ds_in = _coerce_latlon(ds_in)
    try:
        ds_out = rgrd(ds_in, output_chunks=OUT_CHUNKS)
    except TypeError:
        ds_out = rgrd(ds_in).chunk(OUT_CHUNKS)
    key = next((v for v in ds_out.data_vars if {'lat','lon'}.issubset(ds_out[v].dims)), None)
    if key:
        lat_nan = (ds_out[key] == 0).all([d for d in ds_out[key].dims if d != "lat"])
        lon_nan = (ds_out[key] == 0).all([d for d in ds_out[key].dims if d != "lon"])
        ds_out = ds_out.where(~lat_nan & ~lon_nan)
    return ds_out

# ─── 1. Panel ───────────────────────────────────────────────────────────────────
panel = pd.read_stata(PANEL_DTA)
panel["iso"] = collapse_eu(panel["iso"])
panel["year"] = pd.to_numeric(panel["year"], errors="coerce").astype("Int64").dropna().astype(int)
keys = panel[["iso", "adm1_id", "adm2_id", "year"]].drop_duplicates()
spatial_keys = keys[["iso", "adm1_id", "adm2_id"]].drop_duplicates()
years = np.sort(keys["year"].unique())
print(f"[INFO] Panel years {years[0]}–{years[-1]}  ADM2 clusters={len(spatial_keys)}")

# ─── 2. Shapes ──────────────────────────────────────────────────────────────────
print("[INFO] Reading shapefile...")
gdf = gpd.read_file(SHAPEFILE)

# Fix CRS if needed
if gdf.crs is None:
    gdf = gdf.set_crs("EPSG:4326", allow_override=True)
elif gdf.crs.to_epsg() != 4326:
    gdf = gdf.to_crs("EPSG:4326")

# Filter to spatial_keys and ensure consistent dtypes
gdf["iso"] = collapse_eu(gdf["iso"])
gdf["adm1_id"] = gdf["adm1_id"].astype(str)
gdf["adm2_id"] = gdf["adm2_id"].astype(str)
gdf = gdf.merge(spatial_keys, on=["iso", "adm1_id", "adm2_id"], how="inner")

# VALIDITY CHECK: drop null geometries, then apply `make_valid()` only on invalid
gdf = gdf[gdf.geometry.notnull()]
gdf_invalid = gdf[~gdf.is_valid]
if not gdf_invalid.empty:
    print(f"[WARN] Found {len(gdf_invalid)} invalid geometries — fixing them...")
    gdf.loc[gdf_invalid.index, "geometry"] = gdf_invalid.geometry.make_valid()

# FILTER OUT NaN/Inf bounds (Kevin-style fallback)
def has_broken_bounds(geom):
    try:
        b = geom.bounds
        return np.any(np.isnan(b)) or np.any(np.isinf(b))
    except:
        return True

gdf = gdf[~gdf.geometry.apply(has_broken_bounds)].copy()

# Create ADM2 and ADM1 aggregates
gdf_adm2 = gdf.dissolve(by=["iso", "adm1_id", "adm2_id"], as_index=False)
gdf_adm1 = gdf_adm2.dissolve(by=["iso", "adm1_id"], as_index=False)
print(f"[INFO] Shapes aligned: ADM2={len(gdf_adm2)}, ADM1={len(gdf_adm1)}")

# ─── 3. Climate ─────────────────────────────────────────────────────────────────
paths = pd.read_csv(CAR_PATHS_CSV, dtype=str)
row = paths.loc[paths["product"].str.strip() == PRODUCT].iloc[0]
has_precip = (PRODUCT.upper() != "MERRA2" and is_nonempty(row.get("precip_filepath")))

ds_list = []
ds_t = open_any(row["tas_filepath"])
ds_list.append(ds_t[["tas"]])
if has_precip:
    ds_p = open_any(row["precip_filepath"])
    if "prcorr" in ds_p and "pr" not in ds_p:
        ds_p = ds_p.rename({"prcorr": "pr"})
    ds_list.append(ds_p[["pr"]])

ds = xr.merge(ds_list).unify_chunks()
ds = ds.sel(time=ds.time.dt.year.isin(years))
ds["tas"] = (ds["tas"].astype("float32") - 273.15)
if "pr" in ds:
    ds["pr"] = ds["pr"].astype("float32") * 86400.0
ds = ds.chunk({"time": 30, "lat": 120, "lon": 120})

# ─── 4. Regrid + Weights ────────────────────────────────────────────────────────
print("[INFO] Building regridder and computing weightmaps...")
_sample_src = ds.isel(time=0, drop=True)
RGRD = build_regridder(_sample_src)
sample = regrid_xesmf(_sample_src, RGRD).compute()

with xa.set_options(silent=True):
    WM_ADM2 = xa.pixel_overlaps(sample, gdf_adm2)
    WM_ADM1 = xa.pixel_overlaps(sample, gdf_adm1)

# Pick 1 year to debug
yr = 2000
ds_yr = ds.sel(time=ds.time.dt.year == yr)

# Grab tas or pr variable
da = ds_yr["tas"]  # or ds_yr["pr"] if precip

# Just test on the first timestep
da = da.isel(time=0)

# Confirm shape
assert set(da.dims) == {"lat", "lon"}
assert set(WM_ADM2.dims) == {"lat", "lon", "loc"}

agg_test = xr.apply_ufunc(
    numba_aggregate,
    da,
    WM_ADM2,
    input_core_dims=[["lat", "lon"], ["lat", "lon", "loc"]],
    output_core_dims=[["loc"]],
    vectorize=True,
    dask="parallelized",
    output_dtypes=["float32"],
    dask_gufunc_kwargs={"allow_rechunk": True},
)
print(agg_test)