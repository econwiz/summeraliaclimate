#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-
"""
Script A — Per-year ADM2 Aggregation (Carleton-style)

- Reads panel + shapes
- Fixes invalid geometries
- Loads climate (tas + pr)
- Builds weightmap (one-time)
- Loops over panel years:
    * Regrids with cached weights
    * Computes tas/prcp polynomials
    * Aggregates to ADM2
    * Writes one Parquet per year
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import xesmf as xe
import xagg as xa

# ─── CONFIG ────────────────────────────────────────────────────────────────
PRODUCT       = "MERRA2"   # change as needed
BASE          = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data")
CAR_PATHS_CSV = BASE / "car_paths.csv"

PANEL_DTA     = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
SHAPEFILE     = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

OUTDIR_YEARLY = BASE / "obs_yearly"
OUTDIR_YEARLY.mkdir(parents=True, exist_ok=True)

LAT_TARGET = np.arange(-56 + 0.125, 86 - 0.125, 0.25, dtype="float32")
LON_TARGET = np.arange(-180 + 0.125, 180 - 0.125, 0.25, dtype="float32")

# EU collapse
EU_ISO = {
    "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST","FRA","FIN","GBR","GRC",
    "HRV","HUN","IRL","ISL","ITA","LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
    "POL","PRT","ROU","SVK","SVN","SWE","TUR"
}
def collapse_eu(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().where(~s.isin(EU_ISO), "EU")

def tag(name: str) -> str:
    return name.replace("-", "_").upper()

TAG = tag(PRODUCT)

def is_nonempty(x) -> bool:
    return bool(x) and str(x).strip().lower() not in {"", "none", "nan"}

def open_any(path: str) -> xr.Dataset:
    return xr.open_zarr(path, consolidated=False) if path.endswith(".zarr") else xr.open_dataset(path, chunks={"lat":-1,"lon":-1,"time":20})

# ── Regridding helpers ─────────────────────────────────────────────────────
def _coerce_latlon(ds: xr.Dataset) -> xr.Dataset:
    if "latitude" in ds.coords: ds = ds.rename({"latitude":"lat"})
    if "longitude" in ds.coords: ds = ds.rename({"longitude":"lon"})
    if "lon" in ds:
        lon = ds["lon"].values
        if lon.min() >= 0 and lon.max() > 180:
            lon = ((lon+180)%360)-180
            ds = ds.assign_coords(lon=("lon",lon))
        if np.any(np.diff(ds["lon"].values) < 0): ds = ds.sortby("lon")
    if "lat" in ds.coords and np.any(np.diff(ds["lat"].values) < 0): ds = ds.sortby("lat")
    return ds

def make_target_grid() -> xr.Dataset:
    return xr.Dataset(coords={"lat":LAT_TARGET,"lon":LON_TARGET})

def add_bounds_merra(ds_in):
    """Add synthetic lat_b/lon_b cell boundaries for MERRA2 (0.5°x0.625°)."""
    lat = ds_in.lat.values
    lon = ds_in.lon.values
    lat_b = np.concatenate([lat - 0.25, [lat[-1] + 0.25]])
    lon_b = np.concatenate([lon - 0.3125, [lon[-1] + 0.3125]])
    ds_in = ds_in.assign_coords(lat_b=("lat_b", lat_b), lon_b=("lon_b", lon_b))
    return ds_in

def build_regridder(src_like: xr.Dataset, var: str, product: str) -> xe.Regridder:
    src = _coerce_latlon(src_like)
    method = "bilinear"
    if product.upper() == "MERRA2" and var.lower().startswith("pr"):
        src = add_bounds_merra(src)
        method = "conservative"
    dst = make_target_grid()
    # Save weights so we can reuse across runs
    weights_file = BASE / f"weights_{product}_{var}_{method}.nc"
    return xe.Regridder(src, dst, method=method, periodic=False,
                        reuse_weights=os.path.exists(weights_file),
                        filename=str(weights_file))

def regrid_xesmf(ds_in: xr.Dataset, rgrd: xe.Regridder) -> xr.Dataset:
    ds_in = _coerce_latlon(ds_in)
    try:
        return rgrd(ds_in, output_chunks={"time": 20, "lat": 80, "lon": 80})
    except TypeError:
        return rgrd(ds_in).chunk({"time": 20, "lat": 80, "lon": 80})

# ─── 1. Panel ─────────────────────────────────────────────────────────────
panel = pd.read_stata(PANEL_DTA)
panel["iso"] = collapse_eu(panel["iso"])
panel["year"] = pd.to_numeric(panel["year"], errors="coerce").astype("Int64").dropna().astype(int)
keys = panel[["iso","adm1_id","adm2_id","year"]].drop_duplicates()
spatial_keys = keys[["iso","adm1_id","adm2_id"]].drop_duplicates()
years = np.sort(keys["year"].unique())

print(f"[INFO] Panel years {years[0]}–{years[-1]}  ADM2 clusters={len(spatial_keys)}")

# ─── 2. Shapes ─────────────────────────────────────────────────────────────
print("[INFO] Reading shapefile...")
gdf = gpd.read_file(SHAPEFILE).to_crs("EPSG:4326")

# Filter to panel keys
gdf["iso"] = collapse_eu(gdf["iso"])
gdf["adm1_id"] = gdf["adm1_id"].astype(str)
gdf["adm2_id"] = gdf["adm2_id"].astype(str)
gdf = gdf.merge(spatial_keys, on=["iso","adm1_id","adm2_id"], how="inner")

# Fix invalid geometries
gdf = gdf[gdf.geometry.notnull()]
gdf_invalid = gdf[~gdf.is_valid]
if not gdf_invalid.empty:
    print(f"[WARN] Found {len(gdf_invalid)} invalid geometries — fixing...")
    gdf.loc[gdf_invalid.index,"geometry"] = gdf_invalid.geometry.make_valid()

gdf_adm2 = gdf.dissolve(by=["iso","adm1_id","adm2_id"], as_index=False)
print(f"[INFO] Shapes aligned: ADM2={len(gdf_adm2)}")

# ─── 3. Climate ───────────────────────────────────────────────────────────
paths = pd.read_csv(CAR_PATHS_CSV, dtype=str)
row = paths.loc[paths["product"].str.strip() == PRODUCT].iloc[0]
has_precip = (PRODUCT.upper()!="MERRA2" and is_nonempty(row.get("precip_filepath")))

ds_list = []
ds_t = open_any(row["tas_filepath"])
ds_list.append(ds_t[["tas"]])

if has_precip:
    ds_p = open_any(row["precip_filepath"])
    if "prcorr" in ds_p and "pr" not in ds_p:
        ds_p = ds_p.rename({"prcorr":"pr"})
    ds_list.append(ds_p[["pr"]])

ds = xr.merge(ds_list).unify_chunks()
ds = ds.sel(time=ds.time.dt.year.isin(years))
available_years = set(ds.time.dt.year.values)
years = [y for y in years if y in available_years]

ds["tas"] = ds["tas"].astype("float32") - 273.15
if "pr" in ds: ds["pr"] = ds["pr"].astype("float32") * 86400.0

ds = ds.chunk({"lat":-1,"lon":-1,"time":20}).unify_chunks()

print("[DEBUG] ds chunking:", {k:v.chunks for k,v in ds.data_vars.items()})
print(ds.time.min().values, ds.time.max().values)

# ─── 4. Regridder + weightmap (built once) ─────────────────────────────────
print("[INFO] Building weightmaps...")

if PRODUCT.upper() == "MERRA2":
    print("[INFO] Skipping regridding for MERRA2 (native grid)")
    sample = ds.isel(time=0, drop=True)
    if isinstance(sample, xr.DataArray):
        sample = sample.to_dataset(name=sample.name or "tas")
    rgrd_tas = None
    rgrd_pr  = None
else:
    _sample_src = ds["tas"].isel(time=0, drop=True).to_dataset(name="tas")
    rgrd_tas = build_regridder(_sample_src, "tas", PRODUCT)
    rgrd_pr  = build_regridder(ds[["pr"]].isel(time=0), "pr", PRODUCT) if "pr" in ds else None
    sample = regrid_xesmf(ds.isel(time=0, drop=True), rgrd_tas).chunk({"lat": -1, "lon": -1})

with xa.set_options(silent=True):
    WM_ADM2 = xa.pixel_overlaps(sample, gdf_adm2)

# ─── 5. Per-year aggregation ──────────────────────────────────────────────
for yr in years:
    print(f"[INFO] Aggregating {yr}")
    ds_yr = ds.sel(time=ds.time.dt.year == yr)

    # Regrid tas
    tas_regr = regrid_xesmf(ds_yr[["tas"]], rgrd_tas) if rgrd_tas else ds_yr[["tas"]]

    # Regrid pr (if present)
    if "pr" in ds_yr:
        pr_regr = regrid_xesmf(ds_yr[["pr"]], rgrd_pr) if rgrd_pr else ds_yr[["pr"]]
        dsy = xr.merge([tas_regr, pr_regr])
    else:
        dsy = tas_regr

    if yr == years[0]:
        print("[DEBUG] yearly dsy chunks:", {k: v.chunks for k, v in dsy.data_vars.items()})

    # Temperature polynomials
    tas = dsy.tas.astype("float32")
    vars_dict = {
        "T1_sum": tas.sum("time", dtype=np.float32),
        "T2_sum": (tas**2).sum("time", dtype=np.float32),
        "T3_sum": (tas**3).sum("time", dtype=np.float32),
        "T4_sum": (tas**4).sum("time", dtype=np.float32),
    }

    # Precipitation polynomials
    if "pr" in dsy:
        pr = dsy.pr.astype("float32")
        vars_dict.update({
            "PR1_sum": pr.sum("time", dtype=np.float32),
            "PR2_sum": (pr**2).sum("time", dtype=np.float32),
            "PR3_sum": (pr**3).sum("time", dtype=np.float32),
            "PR4_sum": (pr**4).sum("time", dtype=np.float32),
        })

    ds_poly = xr.Dataset(vars_dict)
    with xa.set_options(impl="numba", silent=True):
        agg_da = xa.aggregate(ds_poly, WM_ADM2)

    df = (
        agg_da.to_dataframe()
        .reset_index()
        .rename(columns={
            "T1_sum": f"tavg_poly_1_{TAG}",
            "T2_sum": f"tavg_poly_2_{TAG}",
            "T3_sum": f"tavg_poly_3_{TAG}",
            "T4_sum": f"tavg_poly_4_{TAG}",
            "PR1_sum": f"prcp_poly_1_{TAG}",
            "PR2_sum": f"prcp_poly_2_{TAG}",
            "PR3_sum": f"prcp_poly_3_{TAG}",
            "PR4_sum": f"prcp_poly_4_{TAG}",
        })
        .assign(year=int(yr), product=PRODUCT)
    )

    # Ensure PRCP columns exist
    for k in (1, 2, 3, 4):
        col = f"prcp_poly_{k}_{TAG}"
        if col not in df:
            df[col] = 0.0

    out_path = OUTDIR_YEARLY / f"{PRODUCT.replace('-','_')}_{yr}.parquet"
    df.to_parquet(out_path, index=False)
    print(f"✅ Wrote {out_path}")
