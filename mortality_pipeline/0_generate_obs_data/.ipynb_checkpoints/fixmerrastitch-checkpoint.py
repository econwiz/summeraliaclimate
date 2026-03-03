#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-

"""
Stitch temp ADM2-by-year CSV parts and append ADM1 long-run mean (LR),
without re-running the per-year aggregation.
"""

import os
from pathlib import Path
import re
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import xesmf as xe
import xagg as xa
import dask

# ─── Config ────────────────────────────────────────────────────────────────────
PRODUCT = "ERA5-025"  # change if needed: "GMFD", "ERA5-025", "MERRA2", ...
BASE    = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data")
OUTDIR  = BASE / "obs_csvs"
CAR_PATHS_CSV = BASE / "car_paths.csv"
SHAPEFILE = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

OUT_CSV = OUTDIR / f"{PRODUCT.replace('-', '_')}_by_region_year.csv"

# Target grid used by your other scripts
LAT_TARGET = np.arange(-56 + 0.125, 86 - 0.125, 0.25, dtype="float32")
LON_TARGET = np.arange(-180 + 0.125, 180 - 0.125, 0.25, dtype="float32")

# ─── Dask ──────────────────────────────────────────────────────────────────────
dask.config.set({"array.slicing.split_large_chunks": True})
try:
    from pyproj import datadir
    os.environ["PROJ_LIB"] = datadir.get_data_dir()
except Exception:
    pass

# ─── Helpers ───────────────────────────────────────────────────────────────────
EU_ISO = {
    "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST","FRA","FIN","GBR","GRC",
    "HRV","HUN","IRL","ISL","ITA","LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
    "POL","PRT","ROU","SVK","SVN","SWE","TUR"
}
def collapse_eu(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().where(~s.isin(EU_ISO), "EU")

def tag(name: str) -> str:
    return name.replace("-", "_").upper()

def is_nonempty(x) -> bool:
    return bool(x) and str(x).strip().lower() not in {"", "none", "nan"}

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

def build_regridder(src_like: xr.Dataset, product: str) -> xe.Regridder:
    src = _coerce_latlon(src_like)
    # MERRA2 can choke with periodic=True; disable there.
    periodic = product.upper() not in {"MERRA2"}
    return xe.Regridder(
        src, make_target_grid(),
        method="bilinear",
        periodic=periodic,
        ignore_degenerate=True,
    )

def regrid_xesmf(ds_in: xr.Dataset, rgrd: xe.Regridder) -> xr.Dataset:
    ds_in = _coerce_latlon(ds_in)
    try:
        ds_out = rgrd(ds_in, output_chunks={"time": 12, "lat": 80, "lon": 80})
    except TypeError:
        ds_out = rgrd(ds_in).chunk({"time": 12, "lat": 80, "lon": 80})
    return ds_out

def open_any(path: str) -> xr.Dataset:
    # keep time reasonably small; let spatial follow stored chunks
    return xr.open_zarr(path, consolidated=False) if path.endswith(".zarr") \
        else xr.open_dataset(path, chunks={"time": 20})

# ─── 1) Stitch temp CSVs ───────────────────────────────────────────────────────
def stitch_temp_parts(product: str, outdir: Path) -> pd.DataFrame:
    prefix = f"__tmp_{product}"
    rx = re.compile(rf"^{re.escape(prefix)}.*\.csv$", re.IGNORECASE)
    files = sorted([p for p in outdir.iterdir() if p.is_file() and rx.match(p.name)])
    if not files:
        raise FileNotFoundError(f"No temp CSVs like {prefix}*.csv found in {outdir}")
    parts = []
    for p in files:
        df = pd.read_csv(p)
        # normalize IDs
        for k in ("iso","adm1_id","adm2_id"):
            if k in df.columns:
                df[k] = df[k].astype(str).str.strip()
        if "year" in df.columns:
            df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
        if "iso" in df.columns:
            df["iso"] = collapse_eu(df["iso"])
        parts.append(df)
    stitched = pd.concat(parts, ignore_index=True)
    # de-dup on key keeping the last
    key = ["iso","adm1_id","adm2_id","year"]
    if all(k in stitched.columns for k in key):
        stitched = stitched.drop_duplicates(subset=key, keep="last")
    print(f"[INFO] Stitched {len(files)} temp CSVs → {len(stitched):,} rows")
    return stitched

# ─── 2) Compute ADM1 LR and join ───────────────────────────────────────────────
def compute_lr_and_join(product: str, stitched: pd.DataFrame) -> pd.DataFrame:
    TAG = tag(product)

    # Build ADM1 geometries only for present keys
    keycols = ["iso","adm1_id","adm2_id"]
    assert all(k in stitched.columns for k in keycols), "Key columns missing in stitched CSVs."
    spatial_keys = (
        stitched[keycols].dropna()
        .astype({"adm1_id": str, "adm2_id": str})
        .drop_duplicates()
    )

    gdf = gpd.read_file(SHAPEFILE).to_crs("EPSG:4326")
    gdf["iso"] = collapse_eu(gdf["iso"])
    gdf["adm1_id"] = gdf["adm1_id"].astype(str)
    gdf["adm2_id"] = gdf["adm2_id"].astype(str)
    gdf = gdf.merge(spatial_keys, on=["iso","adm1_id","adm2_id"], how="inner")
    gdf = gdf[gdf.geometry.notnull()]
    bad = gdf[~gdf.is_valid]
    if not bad.empty:
        print(f"[WARN] Fixing {len(bad)} invalid geometries…")
        gdf.loc[bad.index, "geometry"] = bad.geometry.make_valid()
    gdf_adm2 = gdf.dissolve(by=["iso","adm1_id","adm2_id"], as_index=False)
    gdf_adm1 = gdf_adm2.dissolve(by=["iso","adm1_id"], as_index=False)

    # Climate paths & open tas
    paths = pd.read_csv(CAR_PATHS_CSV, dtype=str)
    row = paths.loc[paths["product"].str.strip() == product].iloc[0]
    if not is_nonempty(row.get("tas_filepath")):
        raise RuntimeError("tas_filepath missing in car_paths.csv for product " + product)
    ds_t = open_any(row["tas_filepath"])[["tas"]]
    ds_t["tas"] = (ds_t["tas"].astype("float32") - 273.15)

    # LR window → 1981–2010 if present, else all available
    years_all = np.unique(ds_t.time.dt.year.values)
    years_lr = [y for y in years_all if 1981 <= y <= 2010] or list(years_all)
    ds_lr = ds_t.sel(time=ds_t.time.dt.year.isin(years_lr))

    # Regrid once and compute 2-D LR mean (materialize to avoid huge graphs)
    _sample_src = ds_lr["tas"].isel(time=0, drop=True).to_dataset(name="tas")
    rgrd = build_regridder(_sample_src, product)
    ds_lr_r = regrid_xesmf(ds_lr, rgrd).astype("float32")
    tmean_lr = (
        ds_lr_r.tas
        .groupby("time.year").mean("time")
        .mean("year")
        .transpose("lat","lon")
        .astype("float32")
        .compute()
    )

    # Weights on the regridded grid → aggregate to ADM1
    sample = tmean_lr.to_dataset(name="tas")
    with xa.set_options(silent=True):
        WM_ADM1 = xa.pixel_overlaps(sample, gdf_adm1)

    with xa.set_options(impl="numba", silent=True):
        agg_lr = xa.aggregate(tmean_lr.rename("Tmean"), WM_ADM1)

    df_lr = (
        agg_lr.to_dataframe().reset_index()[["iso","adm1_id","Tmean"]]
        .rename(columns={"Tmean": f"lr_tavg_{TAG}_adm1_avg"})
    )

    # Left-join LR onto stitched table (by iso, adm1)
    out = stitched.merge(df_lr, on=["iso","adm1_id"], how="left")
    return out

# ─── Main ──────────────────────────────────────────────────────────────────────
def main():
    stitched = stitch_temp_parts(PRODUCT, OUTDIR)
    out = compute_lr_and_join(PRODUCT, stitched)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"✅ Wrote {OUT_CSV}  (rows={len(out):,})")

if __name__ == "__main__":
    main()
