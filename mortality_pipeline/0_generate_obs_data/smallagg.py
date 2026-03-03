#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8

import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
import geopandas as gpd
import xagg as xa
from pathlib import Path

# --------- paths ---------
BASE       = Path(__file__).resolve().parent
CAR_PATHS  = BASE / "car_paths.csv"
SHAPEFILE  = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
OUTDIR     = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs")
OUTDIR.mkdir(parents=True, exist_ok=True)

# --------- target grid (Kevin-style 0.25°, centers at +0.125) ---------
LAT_TARGET = np.arange(-56 + 0.125, 86 - 0.124, 0.25, dtype="float32")
LON_TARGET = np.arange(-180 + 0.125, 180 - 0.124, 0.25, dtype="float32")

def _coerce_latlon(ds: xr.Dataset) -> xr.Dataset:
    """rename to lat/lon, wrap 0..360→-180..180, sort increasing"""
    rename = {}
    if "latitude" in ds.coords:  rename["latitude"]  = "lat"
    if "longitude" in ds.coords: rename["longitude"] = "lon"
    if "y" in ds.coords and "lat" not in ds.coords: rename["y"] = "lat"
    if "x" in ds.coords and "lon" not in ds.coords: rename["x"] = "lon"
    if rename:
        ds = ds.rename(rename)
    if "lat" not in ds.coords or "lon" not in ds.coords:
        raise ValueError(f"Missing 'lat'/'lon'. Got coords: {list(ds.coords)}")
    lon = ds["lon"].values
    if lon.min() >= 0 and lon.max() <= 360:
        lon = ((lon + 180) % 360) - 180
        ds = ds.assign_coords(lon=("lon", lon))
    if np.any(np.diff(ds["lat"].values) < 0): ds = ds.sortby("lat")
    if np.any(np.diff(ds["lon"].values) < 0): ds = ds.sortby("lon")
    return ds

def make_target_grid() -> xr.Dataset:
    return xr.Dataset(coords={"lat": (["lat"], LAT_TARGET),
                              "lon": (["lon"], LON_TARGET)})

def build_regridder(src_like: xr.Dataset, product_tag: str) -> xe.Regridder:
    src = _coerce_latlon(src_like)
    tgt = make_target_grid()
    # cache weights per product
    wfile = str(BASE / f"weights_{product_tag}_to025.nc")
    return xe.Regridder(src, tgt, method="bilinear", periodic=True,
                        reuse_weights=True, filename=wfile)

def apply_regrid(ds: xr.Dataset, rg: xe.Regridder) -> xr.Dataset:
    ds = _coerce_latlon(ds)
    for v in ds.data_vars:
        if np.issubdtype(ds[v].dtype, np.floating):
            ds[v] = ds[v].astype("float32")
    out = rg(ds)
    # quick sanity
    try:
        vv = list(out.data_vars)[0]
        print(f"[REGRID] NaN% {vv}: {float(out[vv].isnull().mean())*100:.2f}")
    except Exception:
        pass
    return out

def open_any(path: str) -> xr.Dataset:
    return xr.open_zarr(path, consolidated=False) if str(path).endswith(".zarr") else xr.open_dataset(path)

def product_tag(name: str) -> str:
    return name.replace("-", "_").upper()

# --------- shapes (ADM2 & ADM1) ---------
gdf_adm2 = gpd.read_file(SHAPEFILE, engine="fiona").to_crs("EPSG:4326")
gdf_adm2["geometry"] = gdf_adm2.buffer(0)
gdf_adm2 = gdf_adm2[gdf_adm2.is_valid & ~gdf_adm2.is_empty]
gdf_adm1 = gdf_adm2.dissolve(by=["iso","adm1_id"], as_index=False)

# --------- products from car_paths.csv ---------
paths = pd.read_csv(CAR_PATHS, dtype=str)
paths["product"] = paths["product"].str.strip()
products = list(paths["product"].dropna().unique())

# years
years = range(1981, 2014)
clim_start, clim_end = 1981, 2010

# --------- main loop ---------
for product in products:
    TAG = product_tag(product)
    out_csv = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"
    print(f"\n=== {product} ===")

    # row for this product
    if not (paths["product"] == product).any():
        print(f"  [skip] not in car_paths: {product}")
        continue
    row = paths.loc[paths["product"] == product].iloc[0]

    # open tas
    ds_tas = open_any(row["tas_filepath"]).chunk({"time": 20})
    ds_tas = _coerce_latlon(ds_tas)
    if "tas" not in ds_tas.data_vars:
        print("  [skip] no 'tas' variable")
        continue

    # open precip if present; skip product if missing
    has_pr = ("precip_filepath" in row and isinstance(row["precip_filepath"], str)
              and row["precip_filepath"] and str(row["precip_filepath"]).lower() not in ("nan","none"))
    if not has_pr:
        print("  [skip] no precip filepath (per your note: skip models with no precip, e.g., MERRA)")
        continue

    ds_pr = open_any(row["precip_filepath"]).chunk({"time": 20})
    ds_pr = _coerce_latlon(ds_pr)
    pr_var = next((v for v in ("pr","prcorr","prcp","precip","tp") if v in ds_pr.data_vars), None)
    if pr_var is None:
        print("  [skip] precip variable not found in dataset (skipping)")
        continue
    if pr_var == "prcorr" and "pr" not in ds_pr:
        ds_pr = ds_pr.rename({"prcorr": "pr"})
    elif pr_var in ("prcp","precip","tp") and "pr" not in ds_pr:
        ds_pr = ds_pr.rename({pr_var: "pr"})

    # merge & unit conversions (Kevin-style)
    ds = xr.merge([ds_tas[["tas"]], ds_pr[["pr"]]]).chunk({"time": 20})
    ds["tas"] = (ds["tas"] - 273.15).astype("float32")        # K → °C
    ds["tas"].attrs["units"] = "degC"
    ds["pr"]  = (ds["pr"] * 86400.0).astype("float32")        # kg/m^2/s → mm/day
    ds["pr"].attrs["units"] = "mm/day"

    # build regridder & weights on a regridded slice
    rg = build_regridder(ds.isel(time=0, drop=True), TAG)
    sample = apply_regrid(ds.isel(time=0, drop=True), rg)
    with xa.set_options(silent=True):
        wm_adm2 = xa.pixel_overlaps(sample, gdf_adm2)
        wm_adm1 = xa.pixel_overlaps(sample, gdf_adm1)

    # long-run ADM1 mean T (1981–2010) on regridded grid
    ds_lr   = ds.sel(time=slice(f"{clim_start}-01-01", f"{clim_end}-12-31"))
    ds_lr_r = apply_regrid(ds_lr, rg)
    tmean_lr = ds_lr_r.tas.mean("time").rename("Tmean")
    with xa.set_options(impl="numba", silent=True):
        agg_lr = xa.aggregate(tmean_lr, wm_adm1)
    df_lr = (agg_lr.to_dataframe().reset_index()
             .loc[:, ["iso","adm1_id","Tmean"]]
             .rename(columns={"Tmean": f"lr_tavg_{TAG}_adm1_avg"}))

    # per-year polys: sum of daily powers, then spatial aggregation
    rows_out = []
    for yr in years:
        print(f"  year {yr}")
        dsy   = ds.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))
        dsy_r = apply_regrid(dsy, rg)

        t1 = dsy_r.tas.sum("time").rename("T1_sum")
        t2 = (dsy_r.tas**2).sum("time").rename("T2_sum")
        t3 = (dsy_r.tas**3).sum("time").rename("T3_sum")
        t4 = (dsy_r.tas**4).sum("time").rename("T4_sum")

        p1 = dsy_r.pr.sum("time").rename("PR1_sum")
        p2 = (dsy_r.pr**2).sum("time").rename("PR2_sum")
        p3 = (dsy_r.pr**3).sum("time").rename("PR3_sum")
        p4 = (dsy_r.pr**4).sum("time").rename("PR4_sum")

        dsp = xr.merge([t1,t2,t3,t4,p1,p2,p3,p4])

        with xa.set_options(impl="numba", silent=True):
            agg = xa.aggregate(dsp, wm_adm2)

        df = (agg.to_dataframe().reset_index()
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
              .assign(year=yr, product=product))
        df = df.merge(df_lr, on=["iso","adm1_id"], how="left")
        rows_out.append(df)

    out = pd.concat(rows_out, ignore_index=True)
    out.to_csv(out_csv, index=False)
    print(f"  ✓ wrote {out_csv.name} ({len(out):,} rows)")
