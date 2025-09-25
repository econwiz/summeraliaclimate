#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-

"""
This script generates observational temperature and precipitation data at the ADM2 level
from raw climate model files using the xagg package for spatial aggregation.

Steps:
1. Load raw temperature and precipitation data.
2. Preprocess ADM2-level shapefiles and dissolve them to ADM1 level for computing the lr temperature average.
3. Compute long-run average temperature over the entire sample (1981-2013) at the ADM1 level.
4. Loop over each year and compute:
   - Temperature polynomials (sum of daily tas^p for p=1..4)
   - Annual mean precipitation and its polynomials (prcp^p for p=1..4)
5. Spatially aggregate variables to ADM2 using xagg weight maps.
6. Merge in ADM1-level long-run averages for each year.
7. Write results out as a CSV for merging with the mortality panel later.

Inputs:
- car_paths.csv: Contains paths to temperature and precip files for each product.
- mortality_insample_world.shp: Used for spatial aggregation weights.

Output:
- <product>_by_region_year.csv: Contains ADM2-level temperature and precipitation
  variables for each year, formatted for regression prep.

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

# set products, paths, and chunking for regridding (had to use ChatGPT for help regridding so need to look over that segment because I had trouble following). 
PRODUCT       = "MERRA2" #"MERRA2" #"JRA-3Q" #"GMFD" #ERA5-025
BASE = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data")
CAR_PATHS_CSV = BASE / "car_paths.csv"

PANEL_DTA     = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
SHAPEFILE     = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
OUTDIR        = BASE / "obs_csvs"
OUT_CSV       = OUTDIR / f"{PRODUCT.replace('-', '_')}_by_region_year.csv"
OUTDIR.mkdir(parents=True, exist_ok=True)

#target 0.25 x 0.25 grid for regridding 
LAT_TARGET = np.arange(-56 + 0.125, 86 - 0.125, 0.25, dtype="float32")
LON_TARGET = np.arange(-180 + 0.125, 180 - 0.125, 0.25, dtype="float32")

#chunking for smaller time 
OUT_CHUNKS = {"time": 1, "lat": 40, "lon": 40}

#collpasing the EU isos so it matches the shapefile
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
    return xr.open_zarr(path, consolidated=False) if path.endswith(".zarr") else xr.open_dataset(path, chunks={'lat': -1, 'lon': -1, 'time': 20})

#here's the regridding step... I adapted Kevin's from funcs_preprocessing here to regrid the data correctly. 
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

#This step returns the correct target grid using the inputs
def make_target_grid() -> xr.Dataset:
    return xr.Dataset(coords={"lat": LAT_TARGET, "lon": LON_TARGET})

#Builds the regridder from the climate data 
def build_regridder(src_like: xr.Dataset) -> xe.Regridder:
    src = _coerce_latlon(src_like)
    return xe.Regridder(src, make_target_grid(), method="bilinear", periodic=True, ignore_degenerate=False)

#Regrids and checks for zeros (something Kevin did)
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

#Read the mortality panel and get the right keys so that it is easier to match up later on.  
panel = pd.read_stata(PANEL_DTA)

#collapse panel countries in EU
panel["iso"] = collapse_eu(panel["iso"])
panel["year"] = pd.to_numeric(panel["year"], errors="coerce").astype("Int64").dropna().astype(int)

#Pull only unique adm1/adm2/year combinations
keys = panel[["iso", "adm1_id", "adm2_id", "year"]].drop_duplicates()
spatial_keys = keys[["iso", "adm1_id", "adm2_id"]].drop_duplicates()
years = np.sort(keys["year"].unique())

print(f"[INFO] Panel years {years[0]}–{years[-1]}  ADM2 clusters={len(spatial_keys)}")

# 2. Shapes 
print("[INFO] Reading shapefile...")
gdf = gpd.read_file(SHAPEFILE)
gdf = gdf.to_crs("EPSG:4326")

# Filter shapefile to panel keys only and ensure consistent dtypes
gdf["iso"] = collapse_eu(gdf["iso"])
gdf["adm1_id"] = gdf["adm1_id"].astype(str)
gdf["adm2_id"] = gdf["adm2_id"].astype(str)
gdf = gdf.merge(spatial_keys, on=["iso", "adm1_id", "adm2_id"], how="inner")

#Fix any broken/bad geometries using make_valid. 
gdf = gdf[gdf.geometry.notnull()]
gdf_invalid = gdf[~gdf.is_valid]
if not gdf_invalid.empty:
    print(f"[WARN] Found {len(gdf_invalid)} invalid geometries, fixing them...")
    gdf.loc[gdf_invalid.index, "geometry"] = gdf_invalid.geometry.make_valid()

#Now build two separate versions, one at the adm2 level (for temp and precip), and one at adm1 (for lrtemp)
gdf_adm2 = gdf.dissolve(by=["iso", "adm1_id", "adm2_id"], as_index=False)
gdf_adm1 = gdf_adm2.dissolve(by=["iso", "adm1_id"], as_index=False)
print(f"[INFO] Shapes aligned: ADM2={len(gdf_adm2)}, ADM1={len(gdf_adm1)}")

# 3. Load climate data + agg 
paths = pd.read_csv(CAR_PATHS_CSV, dtype=str)
row = paths.loc[paths["product"].str.strip() == PRODUCT].iloc[0]
has_precip = (PRODUCT.upper() != "MERRA2" and is_nonempty(row.get("precip_filepath")))

#load tas and pr
ds_list = []
ds_t = open_any(row["tas_filepath"])
ds_list.append(ds_t[["tas"]])

if has_precip:
    ds_p = open_any(row["precip_filepath"])
    if "prcorr" in ds_p and "pr" not in ds_p:
        ds_p = ds_p.rename({"prcorr": "pr"})
    ds_list.append(ds_p[["pr"]])

# merge temp and precip and make chunks consistent
ds = xr.merge(ds_list).unify_chunks()

# keep only panel years
ds = ds.sel(time=ds.time.dt.year.isin(years))
available_years = set(ds.time.dt.year.values)
years = [y for y in years if y in available_years]

# convert units as float32
ds["tas"] = (ds["tas"].astype("float32") - 273.15)
if "pr" in ds:
    ds["pr"] = ds["pr"].astype("float32") * 86400.0

#Keep same chunking strategy
ds = ds.chunk({"lat": -1, "lon": -1, "time": 20}).unify_chunks()

print("[DEBUG] ds chunking:", {k: v.chunks for k, v in ds.data_vars.items()})
print(ds.time.min().values, ds.time.max().values)


# 4. Regrid and build weights 

print("[INFO] Building regridder and computing weightmaps...")

#Build the regridder from one slice of time 
_sample_src = ds["tas"].isel(time=0, drop=True)
_sample_src = _sample_src.to_dataset(name="tas")
RGRD = build_regridder(_sample_src)

#now regrid that timestep to get the correct grid sample
sample = regrid_xesmf(ds.isel(time=0, drop=True), RGRD).chunk({"lat": -1, "lon": -1})

#build the overlaps for adm1 and 2 using the regridded sample
with xa.set_options(silent=True):
    WM_ADM2 = xa.pixel_overlaps(sample, gdf_adm2)
    WM_ADM1 = xa.pixel_overlaps(sample, gdf_adm1)

#5. Aggregation 
rows = []

for yr in years:
    print(f"[INFO] Aggregating {yr}")

    # Subset year and check for data presence
    ds_yr = ds.sel(time=ds.time.dt.year == yr)
    if ds_yr.sizes["time"] == 0:
        print(f"[WARN] No data available for {yr}, skipping...")
        continue

    #Rebuild the regridder for each year (this might be bloating the process, but I don't know)
    _sample = ds_yr.isel(time=0, drop=True)
    rgrd = build_regridder(_sample)

    #apply regridding on the dataset and chunk
    dsy = regrid_xesmf(ds_yr, rgrd)
    dsy = dsy.chunk({"lat": -1, "lon": -1, "time": 20}).unify_chunks()
    if yr == years[0]:
        print("[DEBUG] yearly dsy chunks:", {k: v.chunks for k, v in dsy.data_vars.items()})



    #gen temperature polynomials
    tas = dsy.tas.astype("float32")
    vars_dict = {
        "T1_sum": tas.sum("time", dtype=np.float32),
        "T2_sum": (tas ** 2).sum("time", dtype=np.float32),
        "T3_sum": (tas ** 3).sum("time", dtype=np.float32),
        "T4_sum": (tas ** 4).sum("time", dtype=np.float32),
    }

    #gen precipitation polynomials
    if "pr" in dsy:
        pr = dsy.pr.astype("float32")
        vars_dict.update({
            "PR1_sum": pr.sum("time", dtype=np.float32),
            "PR2_sum": (pr ** 2).sum("time", dtype=np.float32),
            "PR3_sum": (pr ** 3).sum("time", dtype=np.float32),
            "PR4_sum": (pr ** 4).sum("time", dtype=np.float32),
        })
    
    #build dataset of all the polys
    ds_poly = xr.Dataset(vars_dict)

    #now agg to adm2
    with xa.set_options(impl="numba", silent=True):
        agg_da = xa.aggregate(ds_poly, WM_ADM2)

    #convert to df, rename columns, add the products
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

#concatenate all rows into one df
df_by_year = pd.concat(rows, ignore_index=True)

#6. ADM1 Long-run Mean. Note that this step was a bit memory heavy, so at one point I had to split it into two separate scripts, but it works for GMFD and JRA

years_lr = [y for y in years if 1981 <= y <= 2010] or list(years)

#regrid again to the target grid
ds_lr_r = regrid_xesmf(ds.sel(time=ds.time.dt.year.isin(years_lr)), RGRD)

#compute the yearly mean then average across years 
tmean_y = ds_lr_r.tas.groupby("time.year").mean("time")
tmean_lr = tmean_y.mean("year").rename("Tmean")

tmean_lr = tmean_lr.transpose("lat", "lon")

#aggregate to the ADM1 level (not adm2)
agg_lr = xa.aggregate(tmean_lr, WM_ADM1)
df_lr = agg_lr.to_dataframe().reset_index()[["iso", "adm1_id", "Tmean"]]
df_lr = df_lr.rename(columns={"Tmean": f"lr_tavg_{TAG}_adm1_avg"})

# 7.Panel merge and write 

#make sure all the ids are consistent 
df_by_year["iso"] = collapse_eu(df_by_year["iso"])
df_by_year["adm1_id"] = df_by_year["adm1_id"].astype(str)
df_by_year["adm2_id"] = df_by_year["adm2_id"].astype(str)

#merge the lr means on the temp poly dataframe
out = df_by_year.merge(df_lr, on=["iso", "adm1_id"], how="left")

#make sure all the adm2 keys are present for all
universe = keys[["iso", "adm1_id", "adm2_id", "year"]].drop_duplicates()
out = universe.merge(out, on=["iso", "adm1_id", "adm2_id", "year"], how="left")

#write out the CSV
out.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV.name} with {len(out):,} rows → {OUTDIR}")

