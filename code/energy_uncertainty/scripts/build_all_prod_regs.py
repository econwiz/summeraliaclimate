#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-
"""
build_PRODUCT_country_climate_multi.py

For each daily climate product listed in car_paths.csv
(ERA5-025, GMFD, JRA-3Q, MERRA2), construct population-weighted
country–year climate variables analogous to the Rode/Carleton TINV
inputs: temp*, precip*, polyAbove*, polyBelow*, hdd20, cdd20.
"""

import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray as rxr
import xagg as xa
import xarray as xr

# ---------------------------------------------------------------------
# Paths / config
# ---------------------------------------------------------------------

YEAR_MIN, YEAR_MAX = 1971, 2010

DATA_DIR = Path(os.environ["DATA"])
WORLD_SHP = DATA_DIR / "shapefiles" / "WORLD" / "gadm28_adm0.shp"

GPW_PATH = Path(
    "/user/ab5405/summeraliaclimate/code/energy_consumption/raw_data/"
    "gpw_adminpoints/gpw_v4_population_count_rev11_2000_15_min.tif"
)

OUT_BASE = Path(
    "/user/ab5405/summeraliaclimate/code/energy_uncertainty/data/country_climate"
)

ISO_COL = "ISO"
TAS_VAR_NAME = "tas"
PR_VAR_NAME = "pr"

SUFFIX_MAP = {
    "ERA5-025": "ERA5",
    "GMFD": "GMFD",
    "JRA-3Q": "JRA_3Q",
    "MERRA2": "MERRA2",
}

#Helpers

def open_any(path: Path) -> xr.Dataset:
    path = Path(path)
    if path.suffix == ".zarr":
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path, chunks={"lat": -1, "lon": -1, "time": 20})


def coerce_lon_lat(ds: xr.Dataset) -> xr.Dataset:
    """Standardize spatial coordinates to lat/lon in [-180, 180], sorted."""
    rename = {}
    if "latitude" in ds.coords:
        rename["latitude"] = "lat"
    if "longitude" in ds.coords:
        rename["longitude"] = "lon"
    if rename:
        ds = ds.rename(rename)

    if "lon" in ds.coords:
        lon = ds["lon"].values
        if lon.min() >= 0 and lon.max() > 180:
            lon = ((lon + 180) % 360) - 180
            ds = ds.assign_coords(lon=("lon", lon))
        if np.any(np.diff(ds["lon"].values) < 0):
            ds = ds.sortby("lon")

    if "lat" in ds.coords and np.any(np.diff(ds["lat"].values) < 0):
        ds = ds.sortby("lat")

    return ds


def compute_daily_transforms(ds: xr.Dataset) -> xr.Dataset:
    """Daily nonlinear transforms for temperature and precipitation."""
    tas_c = ds[TAS_VAR_NAME] - 273.15  #K to C

    T1 = tas_c
    T2 = tas_c ** 2
    T3 = tas_c ** 3
    T4 = tas_c ** 4

    HDD20 = xr.where(tas_c < 20, 20 - tas_c, 0)
    CDD20 = xr.where(tas_c >= 20, tas_c - 20, 0)

    ABV1 = xr.where(tas_c >= 20, tas_c - 20, 0)
    ABV2 = xr.where(tas_c >= 20, tas_c ** 2 - 20 ** 2, 0)
    ABV3 = xr.where(tas_c >= 20, tas_c ** 3 - 20 ** 3, 0)
    ABV4 = xr.where(tas_c >= 20, tas_c ** 4 - 20 ** 4, 0)

    BLW1 = xr.where(tas_c < 20, 20 - tas_c, 0)
    BLW2 = xr.where(tas_c < 20, 20 ** 2 - tas_c ** 2, 0)
    BLW3 = xr.where(tas_c < 20, 20 ** 3 - tas_c ** 3, 0)
    BLW4 = xr.where(tas_c < 20, 20 ** 4 - tas_c ** 4, 0)

    data_vars = dict(
        T1=T1,
        T2=T2,
        T3=T3,
        T4=T4,
        HDD20=HDD20,
        CDD20=CDD20,
        ABV1=ABV1,
        ABV2=ABV2,
        ABV3=ABV3,
        ABV4=ABV4,
        BLW1=BLW1,
        BLW2=BLW2,
        BLW3=BLW3,
        BLW4=BLW4,
    )

    if PR_VAR_NAME is not None and PR_VAR_NAME in ds:
        pr = ds[PR_VAR_NAME]
        data_vars.update(PR1=pr, PR2=pr ** 2)

    return xr.Dataset(data_vars, coords=ds.coords)


def ensure_2d_lat_lon(da: xr.DataArray) -> xr.DataArray:
    """Collapse to a 2D (lat, lon) field, summing over time and averaging other dims."""
    if "time" in da.dims:
        da = da.sum("time", skipna=True)

    for dim in list(da.dims):
        if dim not in ("lat", "lon"):
            da = da.mean(dim)

    if "lat" in da.dims:
        da = da.sortby("lat")
    if "lon" in da.dims:
        da = da.sortby("lon")

    return da


def build_weightmap(sample: xr.DataArray, world: gpd.GeoDataFrame):
    """Compute WORLD × grid weightmap using xagg."""
    world = world.reset_index(drop=True).reset_index(names="shape_id")
    world = world.set_index("shape_id")

    if "time" in sample.dims:
        sample = sample.isel(time=0, drop=True)

    with xa.set_options(silent=True):
        wm = xa.pixel_overlaps(sample, world)

    return wm, world.reset_index()


def load_and_prepare_gpw() -> xr.DataArray:
    """Load GPW GeoTIFF → pop(lat, lon) with proper CRS and fill values cleaned."""
    gpw = rxr.open_rasterio(GPW_PATH)
    gpw = gpw.squeeze("band", drop=True)
    gpw.name = "pop"

    fill = gpw.attrs.get("_FillValue", None)
    if fill is not None:
        gpw = gpw.where(gpw != fill, 0)

    gpw = gpw.rename({"y": "lat", "x": "lon"})
    gpw = gpw.sortby("lat").sortby("lon")

    if gpw.rio.crs is None:
        gpw = gpw.rio.write_crs("EPSG:4326", inplace=True)

    return gpw


def build_product_country_climate(
    product_label: str,
    tas_path: Path,
    precip_path: Path | None,
) -> None:
    """Build country–year climate CSV for a single product."""
    suffix = SUFFIX_MAP.get(product_label, product_label.replace("-", "_"))
    print(f"\n[INFO] Processing product: {product_label} (suffix={suffix})")

    out_dir = OUT_BASE / suffix
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"{suffix}_country_climate_{YEAR_MIN}_{YEAR_MAX}.csv"

    print("[INFO]", product_label, "Loading WORLD shapefile...")
    world = gpd.read_file(WORLD_SHP).to_crs("EPSG:4326")
    if ISO_COL not in world.columns:
        raise ValueError(f"{ISO_COL} not in WORLD columns: {world.columns}")

    print("[INFO]", product_label, "Loading climate...")
    ds_t = open_any(tas_path)
    if precip_path is not None:
        ds_p = open_any(precip_path)
        ds = xr.merge([ds_t, ds_p])
    else:
        ds = ds_t

    ds = coerce_lon_lat(ds)
    ds = ds.sel(time=slice(f"{YEAR_MIN}-01-01", f"{YEAR_MAX}-12-31"))
    years = np.unique(ds.time.dt.year.values)
    print(f"[INFO] {product_label} Years in climate: {years.min()}–{years.max()}")

    print("[INFO]", product_label, "Computing daily transforms...")
    ds_daily = compute_daily_transforms(ds)

    print("[INFO]", product_label, "Loading and reprojecting GPW population...")
    gpw = load_and_prepare_gpw()

    sample = ds_daily["T1"].isel(time=0, drop=True)
    sample = sample.rio.write_crs("EPSG:4326")
    sample = sample.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

    gpw = gpw.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
    gpw_on_grid = gpw.rio.reproject_match(sample).fillna(0).astype("float32")

    if "y" in gpw_on_grid.dims or "x" in gpw_on_grid.dims:
        gpw_on_grid = gpw_on_grid.rename({"y": "lat", "x": "lon"})
    gpw_on_grid = gpw_on_grid.sortby("lat").sortby("lon")

    print("[INFO]", product_label, "Building WORLD × grid weightmap...")
    wm, world_with_id = build_weightmap(sample, world)
    shape_to_iso = world_with_id.set_index("shape_id")[ISO_COL]

    # Country populations
    print("[INFO]", product_label, "Aggregating GPW to country populations...")
    ds_pop = xr.Dataset({"pop": gpw_on_grid}).load()
    with xa.set_options(silent=True):
        agg_pop = xa.aggregate(ds_pop, wm)

    df_pop = agg_pop.to_dataframe().reset_index()
    idx_pop = "shape_id" if "shape_id" in df_pop.columns else df_pop.columns[0]
    df_pop["iso"] = df_pop[idx_pop].map(shape_to_iso)

    pop_by_iso = df_pop.groupby("iso")["pop"].sum()
    pop_by_iso = pop_by_iso[pop_by_iso > 0]

    all_year_rows: list[pd.DataFrame] = []

    # Year loop
    for yr in years:
        print(f"[INFO] {product_label} Aggregating {int(yr)}")
        ds_y = ds_daily.sel(time=ds_daily.time.dt.year == yr)

        ds_sum = xr.Dataset()
        for name, da in ds_y.data_vars.items():
            da_2d = ensure_2d_lat_lon(da)
            if set(da_2d.dims) != {"lat", "lon"}:
                raise RuntimeError(
                    f"{product_label} {name}: still not 2D lat/lon, dims={da_2d.dims}"
                )
            ds_sum[name] = da_2d

        df_year = pop_by_iso.rename("pop_country").to_frame().reset_index()
        df_year["year"] = int(yr)

        for name, da in ds_sum.data_vars.items():
            da2 = da.astype("float32").broadcast_like(gpw_on_grid)
            num = gpw_on_grid * da2
            ds_num = xr.Dataset({"num": num}).load()

            with xa.set_options(silent=True):
                agg_num = xa.aggregate(ds_num, wm)

            df_num = agg_num.to_dataframe().reset_index()
            idx_num = "shape_id" if "shape_id" in df_num.columns else df_num.columns[0]
            df_num["iso"] = df_num[idx_num].map(shape_to_iso)

            num_by_iso = df_num.groupby("iso")["num"].sum()
            df_year[name] = df_year["iso"].map(num_by_iso) / df_year["pop_country"]

        all_year_rows.append(df_year)

    df = pd.concat(all_year_rows, ignore_index=True)

    rename_map = {
        "T1": "temp1_{s}",
        "T2": "temp2_{s}",
        "T3": "temp3_{s}",
        "T4": "temp4_{s}",
        "HDD20": "hdd20_{s}",
        "CDD20": "cdd20_{s}",
        "ABV1": "polyAbove1_{s}",
        "ABV2": "polyAbove2_{s}",
        "ABV3": "polyAbove3_{s}",
        "ABV4": "polyAbove4_{s}",
        "BLW1": "polyBelow1_{s}",
        "BLW2": "polyBelow2_{s}",
        "BLW3": "polyBelow3_{s}",
        "BLW4": "polyBelow4_{s}",
        "PR1": "precip1_{s}",
        "PR2": "precip2_{s}",
    }
    rename_map = {k: v.format(s=suffix) for k, v in rename_map.items()}

    keep_cols = ["iso", "year"] + [k for k in rename_map if k in df.columns]
    df = df[keep_cols].rename(columns=rename_map)
    df = df.sort_values(["iso", "year"]).reset_index(drop=True)

    # Long-run HDD/CDD (time-invariant)
    hdd_col = f"hdd20_{suffix}"
    cdd_col = f"cdd20_{suffix}"

    if hdd_col in df.columns:
        hdd_tinv = df.groupby("iso")[hdd_col].mean()
        df[f"hdd20_TINV_{suffix}"] = df["iso"].map(hdd_tinv)

    if cdd_col in df.columns:
        cdd_tinv = df.groupby("iso")[cdd_col].mean()
        df[f"cdd20_TINV_{suffix}"] = df["iso"].map(cdd_tinv)

    print(f"[INFO] {product_label} Writing {out_csv}")
    df.to_csv(out_csv, index=False)
    print(f"[DONE] {product_label}")

#main driver

def main() -> None:
    car_paths_path = (
        "/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/car_paths.csv"
    )
    print("[INFO] Reading", car_paths_path)
    car_paths = pd.read_csv(car_paths_path)
    car_paths = car_paths.dropna(subset=["product"])

    for _, row in car_paths.iterrows():
        product_label = str(row["product"]).strip()

        if product_label != "MERRA2":
            continue

        tas_path = Path(str(row["tas_filepath"]).strip())
        precip_raw = row.get("precip_filepath", None)
        precip_path = (
            Path(precip_raw.strip())
            if isinstance(precip_raw, str) and precip_raw.strip() != ""
            else None
        )

        build_product_country_climate(
            product_label=product_label,
            tas_path=tas_path,
            precip_path=precip_path,
        )


if __name__ == "__main__":
    main()
