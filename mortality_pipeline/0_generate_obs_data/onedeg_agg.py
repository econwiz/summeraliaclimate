#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-

"""
FLEXIBLE RESOLUTION + MULTI-PRODUCT VERSION

Automatically processes ALL products in a single run.

This script generates observational temperature and precipitation data at the ADM2 level
from raw climate model files using the xagg package for spatial aggregation.

KEY DIFFERENCES FROM ENERGY PIPELINE:
- Uses DAILY precipitation (not monthly!)
- Sums daily precip to get annual totals
- Different specification than Carleton et al. energy paper

MODIFICATIONS:
- Added RESOLUTION parameter (025deg or 1deg)
- Added SKIP_REGRID flag for native grid
- Added PRODUCTS list to process multiple products automatically
- Output directory includes resolution for comparison
- Skip products if output already exists

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

# ============================================================================
# CONFIGURATION - EDIT THESE
# ============================================================================

# Products to process (will loop through all)
PRODUCTS = ["GMFD", "ERA5-025", "JRA-3Q", "MERRA2"]

RESOLUTION    = "1deg"   # "025deg" or "1deg"
SKIP_REGRID   = True     # True if source data already at target resolution

# ============================================================================
# PATHS (shared across all products)
# ============================================================================

BASE = Path("/user/ab5405/summeraliaclimate/code/mortality_pipeline/0_generate_obs_data")
CAR_PATHS_CSV = BASE / "car_paths.csv"

PANEL_DTA     = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
SHAPEFILE     = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

# Output includes resolution tag
OUTDIR_BASE = BASE / "obs_csvs" / RESOLUTION
OUTDIR_BASE.mkdir(parents=True, exist_ok=True)

# ============================================================================
# TARGET GRID (shared)
# ============================================================================

if RESOLUTION == "1deg":
    LAT_TARGET = np.arange(-56 + 0.5, 86 - 0.5, 1.0, dtype="float32")
    LON_TARGET = np.arange(-180 + 0.5, 180 - 0.5, 1.0, dtype="float32")
elif RESOLUTION == "025deg":
    LAT_TARGET = np.arange(-56 + 0.125, 86 - 0.125, 0.25, dtype="float32")
    LON_TARGET = np.arange(-180 + 0.125, 180 - 0.125, 0.25, dtype="float32")
else:
    raise ValueError(f"Unknown RESOLUTION: {RESOLUTION}")

OUT_CHUNKS = {"time": 1, "lat": 40, "lon": 40}

print(f"\n{'='*70}")
print(f"BATCH CONFIGURATION:")
print(f"  Products: {', '.join(PRODUCTS)}")
print(f"  Resolution: {RESOLUTION} ({len(LAT_TARGET)} lats × {len(LON_TARGET)} lons)")
print(f"  Skip Regrid: {SKIP_REGRID}")
print(f"  Output dir: {OUTDIR_BASE}")
print(f"{'='*70}\n")

# ============================================================================
# EU COLLAPSE
# ============================================================================

EU_ISO = {
    "AUT", "BEL", "BGR", "CHE", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FRA", "FIN",
    "GBR", "GRC", "HRV", "HUN", "IRL", "ISL", "ITA", "LIE", "LTU", "LUX", "LVA", "MKD",
    "MLT", "MNE", "NLD", "NOR", "POL", "PRT", "ROU", "SVK", "SVN", "SWE", "TUR"
}

def collapse_eu(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().where(~s.isin(EU_ISO), "EU")

def is_nonempty(x) -> bool:
    return bool(x) and str(x).strip().lower() not in {"", "none", "nan"}

def open_any(path: str) -> xr.Dataset:
    return xr.open_zarr(path, consolidated=False) if path.endswith(".zarr") else xr.open_dataset(path, chunks={'lat': -1, 'lon': -1, 'time': 20})

# ============================================================================
# REGRIDDING FUNCTIONS
# ============================================================================

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

# ============================================================================
# LOAD PANEL AND SHAPEFILE (once, shared across products)
# ============================================================================

print("[INFO] Loading mortality panel...")
panel = pd.read_stata(PANEL_DTA)
panel["iso"] = collapse_eu(panel["iso"])
panel["year"] = pd.to_numeric(panel["year"], errors="coerce").astype("Int64").dropna().astype(int)

keys = panel[["iso", "adm1_id", "adm2_id", "year"]].drop_duplicates()
spatial_keys = keys[["iso", "adm1_id", "adm2_id"]].drop_duplicates()
years = np.sort(keys["year"].unique())

print(f"[INFO] Panel years {years[0]}–{years[-1]}  ADM2 clusters={len(spatial_keys)}")

print("[INFO] Loading shapefile...")
gdf = gpd.read_file(SHAPEFILE)
gdf = gdf.to_crs("EPSG:4326")

gdf["iso"] = collapse_eu(gdf["iso"])
gdf["adm1_id"] = gdf["adm1_id"].astype(str)
gdf["adm2_id"] = gdf["adm2_id"].astype(str)
gdf = gdf.merge(spatial_keys, on=["iso", "adm1_id", "adm2_id"], how="inner")

gdf = gdf[gdf.geometry.notnull()]
gdf_invalid = gdf[~gdf.is_valid]
if not gdf_invalid.empty:
    print(f"[WARN] Found {len(gdf_invalid)} invalid geometries, fixing...")
    gdf.loc[gdf_invalid.index, "geometry"] = gdf_invalid.geometry.make_valid()

gdf_adm2 = gdf.dissolve(by=["iso", "adm1_id", "adm2_id"], as_index=False)
gdf_adm1 = gdf_adm2.dissolve(by=["iso", "adm1_id"], as_index=False)
print(f"[INFO] Shapes: ADM2={len(gdf_adm2)}, ADM1={len(gdf_adm1)}\n")

# ============================================================================
# MAIN PRODUCT LOOP
# ============================================================================

for PRODUCT in PRODUCTS:
    
    print(f"\n{'#'*70}")
    print(f"# PROCESSING PRODUCT: {PRODUCT}")
    print(f"{'#'*70}\n")
    
    TAG = PRODUCT.replace("-", "_").upper()
    OUT_CSV = OUTDIR_BASE / f"{TAG}_by_region_year.csv"
    
    # Skip if output already exists
    if OUT_CSV.exists():
        print(f"[SKIP] Output already exists: {OUT_CSV}\n")
        continue
    
    try:
        # ====================================================================
        # LOAD CLIMATE DATA FOR THIS PRODUCT
        # ====================================================================
        
        print(f"[INFO] Loading climate data for {PRODUCT}...")
        paths = pd.read_csv(CAR_PATHS_CSV, dtype=str)
        row = paths.loc[paths["product"].str.strip() == PRODUCT].iloc[0]
        has_precip = is_nonempty(row.get("precip_filepath"))

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
        available_years = set(ds.time.dt.year.values)
        years_product = [y for y in years if y in available_years]

        # Convert units
        ds["tas"] = (ds["tas"].astype("float32") - 273.15)
        if "pr" in ds:
            ds["pr"] = ds["pr"].astype("float32") * 86400.0

        ds = ds.chunk({"lat": -1, "lon": -1, "time": 20}).unify_chunks()

        print(f"[INFO] Time range: {ds.time.min().values} to {ds.time.max().values}")

        # ====================================================================
        # BUILD REGRIDDER/WEIGHTS (or skip if native grid)
        # ====================================================================

        if SKIP_REGRID:
            print(f"[INFO] SKIP_REGRID=True, using native grid")
            sample = ds.isel(time=0, drop=True).chunk({"lat": -1, "lon": -1})
            RGRD = None
            
            lat_res = np.diff(ds.lat.values).mean()
            lon_res = np.diff(ds.lon.values).mean()
            print(f"[INFO] Source grid: {lat_res:.4f}° lat × {lon_res:.4f}° lon")
            
            if RESOLUTION == "1deg" and not (0.9 < lat_res < 1.1):
                print(f"[WARN] Expected 1° but got {lat_res:.4f}° - set SKIP_REGRID=False!")
            elif RESOLUTION == "025deg" and not (0.2 < lat_res < 0.3):
                print(f"[WARN] Expected 0.25° but got {lat_res:.4f}° - set SKIP_REGRID=False!")
        else:
            print(f"[INFO] Building regridder to {RESOLUTION}...")
            _sample_src = ds["tas"].isel(time=0, drop=True).to_dataset(name="tas")
            RGRD = build_regridder(_sample_src)
            sample = regrid_xesmf(ds.isel(time=0, drop=True), RGRD).chunk({"lat": -1, "lon": -1})

        with xa.set_options(silent=True):
            WM_ADM2 = xa.pixel_overlaps(sample, gdf_adm2)
            WM_ADM1 = xa.pixel_overlaps(sample, gdf_adm1)

        print("[INFO] Weight maps built")

        # ====================================================================
        # YEARLY AGGREGATION
        # ====================================================================

        rows = []

        for yr in years_product:
            print(f"[INFO] Aggregating {yr}")

            ds_yr = ds.sel(time=ds.time.dt.year == yr)
            if ds_yr.sizes["time"] == 0:
                print(f"[WARN] No data for {yr}, skipping...")
                continue

            # Apply regridding if needed
            if SKIP_REGRID:
                dsy = ds_yr.chunk({"lat": -1, "lon": -1, "time": 20}).unify_chunks()
            else:
                _sample = ds_yr.isel(time=0, drop=True)
                rgrd = build_regridder(_sample)
                dsy = regrid_xesmf(ds_yr, rgrd)
                dsy = dsy.chunk({"lat": -1, "lon": -1, "time": 20}).unify_chunks()

            # Temperature polynomials
            tas = dsy.tas.astype("float32")
            vars_dict = {
                "T1_sum": tas.sum("time", dtype=np.float32),
                "T2_sum": (tas ** 2).sum("time", dtype=np.float32),
                "T3_sum": (tas ** 3).sum("time", dtype=np.float32),
                "T4_sum": (tas ** 4).sum("time", dtype=np.float32),
            }

            # Precipitation polynomials (DAILY!)
            if "pr" in dsy:
                pr = dsy.pr.astype("float32")
                vars_dict.update({
                    "PR1_sum": pr.sum("time", dtype=np.float32),
                    "PR2_sum": (pr ** 2).sum("time", dtype=np.float32),
                    "PR3_sum": (pr ** 3).sum("time", dtype=np.float32),
                    "PR4_sum": (pr ** 4).sum("time", dtype=np.float32),
                })
            
            ds_poly = xr.Dataset(vars_dict)

            with xa.set_options(impl="numba", silent=True):
                agg_da = xa.aggregate(ds_poly, WM_ADM2)

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

            rows.append(df)

        df_by_year = pd.concat(rows, ignore_index=True)

        # ====================================================================
        # ADM1 LONG-RUN MEAN
        # ====================================================================

        print("[INFO] Computing ADM1 long-run temperature...")

        years_lr = [y for y in years_product if 1981 <= y <= 2010] or list(years_product)

        if SKIP_REGRID:
            ds_lr_r = ds.sel(time=ds.time.dt.year.isin(years_lr))
        else:
            ds_lr_r = regrid_xesmf(ds.sel(time=ds.time.dt.year.isin(years_lr)), RGRD)

        tmean_y = ds_lr_r.tas.groupby("time.year").mean("time")
        tmean_lr = tmean_y.mean("year").rename("Tmean")
        tmean_lr = tmean_lr.transpose("lat", "lon")

        agg_lr = xa.aggregate(tmean_lr, WM_ADM1)
        df_lr = agg_lr.to_dataframe().reset_index()[["iso", "adm1_id", "Tmean"]]
        df_lr = df_lr.rename(columns={"Tmean": f"lr_tavg_{TAG}_adm1_avg"})

        # ====================================================================
        # MERGE AND WRITE
        # ====================================================================

        print("[INFO] Merging and writing output...")

        df_by_year["iso"] = collapse_eu(df_by_year["iso"])
        df_by_year["adm1_id"] = df_by_year["adm1_id"].astype(str)
        df_by_year["adm2_id"] = df_by_year["adm2_id"].astype(str)

        out = df_by_year.merge(df_lr, on=["iso", "adm1_id"], how="left")

        universe = keys[["iso", "adm1_id", "adm2_id", "year"]].drop_duplicates()
        out = universe.merge(out, on=["iso", "adm1_id", "adm2_id", "year"], how="left")

        out.to_csv(OUT_CSV, index=False)

        print(f"\n{'='*70}")
        print(f"✅ SUCCESS - {PRODUCT}")
        print(f"Wrote {OUT_CSV.name}")
        print(f"  Rows: {len(out):,}")
        print(f"{'='*70}\n")
        
    except Exception as e:
        print(f"\n{'!'*70}")
        print(f"❌ ERROR processing {PRODUCT}:")
        print(f"  {str(e)}")
        print(f"{'!'*70}\n")
        continue

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n{'='*70}")
print(f"BATCH PROCESSING COMPLETE")
print(f"Resolution: {RESOLUTION}")
print(f"Output directory: {OUTDIR_BASE}")
print(f"\nGenerated files:")
for prod in PRODUCTS:
    tag = prod.replace("-", "_").upper()
    csv_file = OUTDIR_BASE / f"{tag}_by_region_year.csv"
    status = "✅" if csv_file.exists() else "❌"
    print(f"  {status} {csv_file.name}")
print(f"{'='*70}\n")