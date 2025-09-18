#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-

"""
Climate Aggregation Script — ADM2-by-Year + ADM1 Long-Run Mean (Carleton-style)

- Handles GMFD, ERA5-025, etc.
- Aggregates tas and precip (power 1–4) over ADM2-by-year
- Computes ADM1 long-run mean temperature
- Panel-anchored output with EU collapse
- Memory-safe: single regridder; per-year streaming to CSV; restricted merges
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import xesmf as xe
import xagg as xa
import dask
import gc

# ─── Dask setup ────────────────────────────────────────────────────────────────
dask.config.set({"array.slicing.split_large_chunks": True})
try:
    from pyproj import datadir
    os.environ["PROJ_LIB"] = datadir.get_data_dir()
except Exception:
    pass

# ─── User Settings ─────────────────────────────────────────────────────────────
PRODUCT       = "GMFD"  # "MERRA2" / "JRA-3Q" / "ERA5-025" etc.
BASE = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data")
CAR_PATHS_CSV = BASE / "car_paths.csv"

PANEL_DTA     = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
SHAPEFILE     = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
OUTDIR        = BASE / "obs_csvs"
OUT_CSV       = OUTDIR / f"{PRODUCT.replace('-', '_')}_by_region_year.csv"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Optional: cache xESMF weights to disk so re-runs are instant
WEIGHTS_FILE  = OUTDIR / f"weights_{PRODUCT.replace('-', '_')}_0p25.nc"

LAT_TARGET = np.arange(-56 + 0.125, 86 - 0.125, 0.25, dtype="float32")
LON_TARGET = np.arange(-180 + 0.125, 180 - 0.125, 0.25, dtype="float32")
OUT_CHUNKS = {"lat": -1, "lon": -1, "time": 20}

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
    # unifying initial chunking keeps graphs small
    if str(path).endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path, chunks={"lat": -1, "lon": -1, "time": 20})

# ─── xESMF helpers ─────────────────────────────────────────────────────────────
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

def build_regridder(src_like: xr.Dataset, weights_file: Path | None = None) -> xe.Regridder:
    src = _coerce_latlon(src_like)
    kwargs = dict(method="bilinear", periodic=True, ignore_degenerate=False)
    if weights_file is not None:
        kwargs.update({
            "filename": str(weights_file),
            "reuse_weights": weights_file.exists()
        })
    return xe.Regridder(src, make_target_grid(), **kwargs)

def regrid_xesmf(ds_in: xr.Dataset, rgrd: xe.Regridder) -> xr.Dataset:
    ds_in = _coerce_latlon(ds_in)
    try:
        ds_out = rgrd(ds_in, output_chunks=OUT_CHUNKS)
    except TypeError:
        ds_out = rgrd(ds_in).chunk(OUT_CHUNKS)
    # mask degenerate rows/cols if any (rare but helps)
    key = next((v for v in ds_out.data_vars if {"lat","lon"}.issubset(ds_out[v].dims)), None)
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
years_all = np.sort(keys["year"].unique())
print(f"[INFO] Panel years {years_all[0]}–{years_all[-1]}  ADM2 clusters={len(spatial_keys)}")

# ─── 2. Shapes ──────────────────────────────────────────────────────────────────
print("[INFO] Reading shapefile...")
gdf = gpd.read_file(SHAPEFILE).to_crs("EPSG:4326")

# Filter to spatial_keys and ensure consistent dtypes
gdf["iso"] = collapse_eu(gdf["iso"])
gdf["adm1_id"] = gdf["adm1_id"].astype(str)
gdf["adm2_id"] = gdf["adm2_id"].astype(str)
gdf = gdf.merge(spatial_keys, on=["iso", "adm1_id", "adm2_id"], how="inner")

# Validity cleanup
gdf = gdf[gdf.geometry.notnull()]
gdf_invalid = gdf[~gdf.is_valid]
if not gdf_invalid.empty:
    print(f"[WARN] Found {len(gdf_invalid)} invalid geometries — fixing them...")
    gdf.loc[gdf_invalid.index, "geometry"] = gdf_invalid.geometry.make_valid()

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

# Merge + clip to panel years available in the data
ds = xr.merge(ds_list).unify_chunks()
ds = ds.chunk({"lat": -1, "lon": -1, "time": 20}).unify_chunks()

#Keep only panel years; also filter our `years` list to those actually present
ds_years = np.unique(ds.time.dt.year.values)
years = [int(y) for y in years_all if y in ds_years]
ds = ds.sel(time=ds.time.dt.year.isin(years))

# Unit conversions as float32
ds["tas"] = (ds["tas"].astype("float32") - 273.15)
if "pr" in ds:
    ds["pr"] = ds["pr"].astype("float32") * 86400.0  # from kg m-2 s-1 to mm/day

print("[DEBUG] ds chunking:", {k: v.chunks for k, v in ds.data_vars.items()})
print(ds.time.min().values, ds.time.max().values)

#Build a regridder and compute weightmaps once on a representative slice
print("[INFO] Building regridder and computing weightmaps...")
_sample_src = ds["tas"].isel(time=0, drop=True).to_dataset(name="tas")
RGRD = build_regridder(_sample_src, WEIGHTS_FILE)

# Single regrid of one time-slice to establish target grid for weights
sample = regrid_xesmf(ds.isel(time=0, drop=True), RGRD).chunk({"lat": -1, "lon": -1})

with xa.set_options(silent=True):
    WM_ADM2 = xa.pixel_overlaps(sample, gdf_adm2)
    WM_ADM1 = xa.pixel_overlaps(sample, gdf_adm1)

# Annual ADM2 aggregation: polynomial **sums** over days, then spatial average
tmp_csv = OUTDIR / f"__tmp_{TAG}_by_region_year.csv"
if tmp_csv.exists():
    tmp_csv.unlink()
first_write = True

for yr in years:
    print(f"[INFO] Aggregating {yr}")

    ds_yr = ds.sel(time=ds.time.dt.year == yr)
    if ds_yr.sizes.get("time", 0) == 0:
        print(f"[WARN] No data available for {yr}, skipping...")
        continue

    # Regrid with the single, cached regridder
    dsy = regrid_xesmf(ds_yr, RGRD).chunk({"lat": -1, "lon": -1, "time": 20}).unify_chunks()
    if yr == years[0]:
        print("[DEBUG] yearly dsy chunks:", {k: v.chunks for k, v in dsy.data_vars.items()})

    # Temperature polynomials
    tas = dsy.tas.astype("float32")
    vars_dict = {
        "T1_sum": tas.sum("time", dtype=np.float32),
        "T2_sum": (tas ** 2).sum("time", dtype=np.float32),
        "T3_sum": (tas ** 3).sum("time", dtype=np.float32),
        "T4_sum": (tas ** 4).sum("time", dtype=np.float32),
    }

    # Precipitation polynomials
    if "pr" in dsy:
        pr = dsy.pr.astype("float32")
        vars_dict.update({
            "PR1_sum": pr.sum("time", dtype=np.float32),
            "PR2_sum": (pr ** 2).sum("time", dtype=np.float32),
            "PR3_sum": (pr ** 3).sum("time", dtype=np.float32),
            "PR4_sum": (pr ** 4).sum("time", dtype=np.float32),
        })

    ds_poly = xr.Dataset(vars_dict)

    # Spatial aggregation to ADM2
    with xa.set_options(impl="numba", silent=True):
        agg_da = xa.aggregate(ds_poly, WM_ADM2)

    # To DataFrame and stream-append
    df = (
        agg_da.to_dataframe()
        .reset_index()
        .rename(columns={
            "T1_sum":  f"tavg_poly_1_{TAG}",
            "T2_sum":  f"tavg_poly_2_{TAG}",
            "T3_sum":  f"tavg_poly_3_{TAG}",
            "T4_sum":  f"tavg_poly_4_{TAG}",
            "PR1_sum": f"prcp_poly_1_{TAG}",
            "PR2_sum": f"prcp_poly_2_{TAG}",
            "PR3_sum": f"prcp_poly_3_{TAG}",
            "PR4_sum": f"prcp_poly_4_{TAG}",
        })
        .assign(year=int(yr), product=PRODUCT)
    )

    # Ensure PRCP columns exist if precip missing for product
    for k in (1, 2, 3, 4):
        col = f"prcp_poly_{k}_{TAG}"
        if col not in df:
            df[col] = np.float32(0.0)

    # Normalize types + write
    df["iso"] = collapse_eu(df["iso"])
    df["adm1_id"] = df["adm1_id"].astype(str)
    df["adm2_id"] = df["adm2_id"].astype(str)
    num_cols = [c for c in df.columns if c.startswith(("tavg_poly_", "prcp_poly_"))]
    df[num_cols] = df[num_cols].astype("float32")
    df.to_csv(tmp_csv, mode="a", header=first_write, index=False)
    first_write = False

    # Free memory this iteration
    del ds_yr, dsy, tas, ds_poly, agg_da, df
    try:
        del pr
    except NameError:
        pass
    gc.collect()

# Gather the streamed results
df_by_year = pd.read_csv(tmp_csv, low_memory=False)

years_lr = [y for y in years if 1981 <= y <= 2010] or list(years)
ds_lr_r = regrid_xesmf(ds.sel(time=ds.time.dt.year.isin(years_lr)), RGRD)
tmean_y = ds_lr_r.tas.groupby("time.year").mean("time")
tmean_lr = tmean_y.mean("year").rename("Tmean").transpose("lat", "lon")

with xa.set_options(impl="numba", silent=True):
    agg_lr = xa.aggregate(tmean_lr, WM_ADM1)

df_lr = agg_lr.to_dataframe().reset_index()[["iso", "adm1_id", "Tmean"]]
df_lr["iso"] = collapse_eu(df_lr["iso"])
df_lr["adm1_id"] = df_lr["adm1_id"].astype(str)
df_lr = df_lr.rename(columns={"Tmean": f"lr_tavg_{TAG}_adm1_avg"})

df_by_year["iso"] = collapse_eu(df_by_year["iso"])
df_by_year["adm1_id"] = df_by_year["adm1_id"].astype(str)
df_by_year["adm2_id"] = df_by_year["adm2_id"].astype(str)

universe = keys[keys["year"].isin(years)][["iso", "adm1_id", "adm2_id", "year"]].drop_duplicates()

out = df_by_year.merge(df_lr, on=["iso", "adm1_id"], how="left")
out = universe.merge(out, on=["iso", "adm1_id", "adm2_id", "year"], how="left")

out.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV.name} with {len(out):,} rows → {OUTDIR}")

try:
    tmp_csv.unlink()
except Exception:
    pass
