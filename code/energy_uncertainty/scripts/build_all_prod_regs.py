#!/user/ab5405/.conda/envs/hle_iv/bin/python
# -*- coding: utf-8 -*-
"""
build_PRODUCT_country_climate_multi.py

For each climate product listed in car_paths.csv
(ERA5-025, GMFD, JRA-3Q, MERRA2), construct population-weighted
country–year climate variables analogous to the Rode/Carleton TINV
inputs: temp*, precip*, polyAbove*, polyBelow*, hdd20, cdd20.

- Temperature: uses DAILY tas (tas_day_*), builds daily transforms,
  then sums over days within the year.
- Precipitation: uses DAILY precip flux (pr_*), converts
  kg m-2 s-1 -> mm/day, then sums over days within the year.

Inputs:
- car_paths.csv: contains paths to tas and precip files for each product.
- gadm28_adm0.shp: country-level shapefile for spatial aggregation.
- gpw_v4_population_count_rev11_2000_15_min.tif: GPW year-2000 population
  raster used as static aggregation weights (consistent with Rode et al.).

Output:
- <SUFFIX>/<SUFFIX>_country_year_<YEAR_MIN>_<YEAR_MAX>.csv:
  Country-year panel of population-weighted climate variables,
  including temperature polynomials, HDD/CDD, polyAbove/polyBelow,
  their TINV long-run means, and pixel-level CDD/HDD cross terms.
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

# GPW v4 year-2000 population raster used as static aggregation weights.
# Year 2000 is used for all years following the Rode et al. convention
# of fixed population weights to avoid confounding climate aggregation
# with demographic change.
GPW_PATH = Path(
    "/user/ab5405/summeraliaclimate/code/energy_consumption/raw_data/"
    "gpw_adminpoints/gpw_v4_population_count_rev11_2000_15_min.tif"
)

OUT_BASE = Path(
    "/user/ab5405/summeraliaclimate/code/energy_uncertainty/data/country_climate"
)

ISO_COL = "ISO"
TAS_VAR_NAME = "tas"
PRDAY_VAR_NAME = "pr"

# Maps product labels (which may contain hyphens) to the filename/column suffix.
SUFFIX_MAP = {
    "ERA5-025": "ERA5",
    "GMFD": "GMFD",
    "JRA-3Q": "JRA_3Q",
    "MERRA2": "MERRA2",
}

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def open_any(path: Path) -> xr.Dataset:
    path = Path(path)
    if path.suffix == ".zarr":
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path, chunks={"lat": -1, "lon": -1, "time": 20})


def coerce_lon_lat(ds: xr.Dataset) -> xr.Dataset:
    """Standardize coordinate names and ensure longitudes are in [-180, 180]
    and sorted ascending, as required by xESMF and xagg."""
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
    """
    Compute daily nonlinear temperature transforms at each grid cell.

    All transforms use 20°C as the threshold, following Rode et al. and the
    CIL energy demand specification. ABV variables capture exposure above 20°C
    (cooling demand analog) and BLW variables capture exposure below 20°C
    (heating demand analog). HDD20/CDD20 are standard heating/cooling degree days.

    For this pipeline we use daily tas only. Precipitation is handled
    separately from daily pr files, so PRDAY_VAR_NAME is not expected in ds.
    """
    tas_c = ds[TAS_VAR_NAME] - 273.15  # Convert Kelvin to Celsius

    T1 = tas_c
    T2 = tas_c ** 2
    T3 = tas_c ** 3
    T4 = tas_c ** 4

    HDD20 = xr.where(tas_c < 20, 20 - tas_c, 0)
    CDD20 = xr.where(tas_c >= 20, tas_c - 20, 0)

    # polyAbove: polynomial exposure above 20°C, normalized to zero at threshold
    ABV1 = xr.where(tas_c >= 20, tas_c - 20, 0)
    ABV2 = xr.where(tas_c >= 20, tas_c ** 2 - 20 ** 2, 0)
    ABV3 = xr.where(tas_c >= 20, tas_c ** 3 - 20 ** 3, 0)
    ABV4 = xr.where(tas_c >= 20, tas_c ** 4 - 20 ** 4, 0)

    # polyBelow: polynomial exposure below 20°C, normalized to zero at threshold
    BLW1 = xr.where(tas_c < 20, 20 - tas_c, 0)
    BLW2 = xr.where(tas_c < 20, 20 ** 2 - tas_c ** 2, 0)
    BLW3 = xr.where(tas_c < 20, 20 ** 3 - tas_c ** 3, 0)
    BLW4 = xr.where(tas_c < 20, 20 ** 4 - tas_c ** 4, 0)

    data_vars = dict(
        T1=T1, T2=T2, T3=T3, T4=T4,
        HDD20=HDD20, CDD20=CDD20,
        ABV1=ABV1, ABV2=ABV2, ABV3=ABV3, ABV4=ABV4,
        BLW1=BLW1, BLW2=BLW2, BLW3=BLW3, BLW4=BLW4,
    )

    return xr.Dataset(data_vars, coords=ds.coords)


def ensure_2d_lat_lon(da: xr.DataArray, *, time_agg: str | None = "sum") -> xr.DataArray:
    """
    Collapse a DataArray to a 2D (lat, lon) field.

    - If time_agg == "sum": sum over the time dimension.
    - If time_agg == "mean": mean over the time dimension.
    - If time_agg is None: time dimension already collapsed by caller; leave alone.
    Any remaining non-(lat, lon) dims are averaged out.
    """
    if "time" in da.dims and time_agg is not None:
        if time_agg == "sum":
            da = da.sum("time", skipna=True)
        elif time_agg == "mean":
            da = da.mean("time", skipna=True)
        else:
            raise ValueError(f"Unknown time_agg='{time_agg}' (expected 'sum', 'mean', or None).")

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


# ---------------------------------------------------------------------
# Main worker
# ---------------------------------------------------------------------

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
    out_csv = out_dir / f"{suffix}_country_year_{YEAR_MIN}_{YEAR_MAX}.csv"

    print("[INFO]", product_label, "Loading WORLD shapefile")
    world = gpd.read_file(WORLD_SHP).to_crs("EPSG:4326")
    if ISO_COL not in world.columns:
        raise ValueError(f"{ISO_COL} not in WORLD columns: {world.columns}")

    # Load temperature
    print("[INFO]", product_label, "Loading DAILY temperature")
    ds_t = open_any(tas_path)
    ds_t = coerce_lon_lat(ds_t)
    ds_t = ds_t.sel(time=slice(f"{YEAR_MIN}-01-01", f"{YEAR_MAX}-12-31"))

    # Load precipitation
    if precip_path is not None:
        print("[INFO]", product_label, "Loading DAILY precipitation")
        ds_p = open_any(precip_path)
        ds_p = coerce_lon_lat(ds_p)
        ds_p = ds_p.sel(time=slice(f"{YEAR_MIN}-01-01", f"{YEAR_MAX}-12-31"))
    else:
        ds_p = None

    # Determine overlapping years between temperature and precipitation
    years_t = np.unique(ds_t.time.dt.year.values)
    if ds_p is not None:
        years_p = np.unique(ds_p.time.dt.year.values)
        years = np.intersect1d(years_t, years_p)
    else:
        years = years_t

    print(f"[INFO] {product_label} Years in temperature: {years_t.min()}–{years_t.max()}")
    if ds_p is not None:
        print(f"[INFO] {product_label} Years in precipitation: {years_p.min()}–{years_p.max()}")
    print(f"[INFO] {product_label} Using overlapping years: {years.min()}–{years.max()}")

    print("[INFO]", product_label, "Computing daily transforms (temperature)")
    ds_daily = compute_daily_transforms(ds_t)

    print("[INFO]", product_label, "Loading and reprojecting GPW population")
    gpw = load_and_prepare_gpw()

    sample = ds_daily["T1"].isel(time=0, drop=True)
    sample = sample.rio.write_crs("EPSG:4326")
    sample = sample.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

    gpw = gpw.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
    gpw_on_grid = gpw.rio.reproject_match(sample).fillna(0).astype("float32")

    if "y" in gpw_on_grid.dims or "x" in gpw_on_grid.dims:
        gpw_on_grid = gpw_on_grid.rename({"y": "lat", "x": "lon"})
    gpw_on_grid = gpw_on_grid.sortby("lat").sortby("lon")

    print("[INFO]", product_label, "Building WORLD × grid weightmap")
    wm, world_with_id = build_weightmap(sample, world)
    shape_to_iso = world_with_id.set_index("shape_id")[ISO_COL]

    # Aggregate GPW population to country level for use as weights denominator
    print("[INFO]", product_label, "Aggregating GPW to country populations")
    ds_pop = xr.Dataset({"pop": gpw_on_grid}).load()
    with xa.set_options(silent=True):
        agg_pop = xa.aggregate(ds_pop, wm)

    df_pop = agg_pop.to_dataframe().reset_index()
    idx_pop = "shape_id" if "shape_id" in df_pop.columns else df_pop.columns[0]
    df_pop["iso"] = df_pop[idx_pop].map(shape_to_iso)

    pop_by_iso = df_pop.groupby("iso")["pop"].sum()
    pop_by_iso = pop_by_iso[pop_by_iso > 0]

    # -------------------------------------------------------------------------
    # Compute pixel-level TINV (time-invariant long-run mean HDD/CDD).
    # Done once here in a single vectorized groupby before the year loop,
    # rather than re-iterating over years. Used as the heterogeneity moderator
    # in the CIL cross-term specification: polyAbove × CDD_TINV, polyBelow × HDD_TINV.
    # Cross products must be formed at the pixel level before aggregation —
    # computing TINV at the country level first would give a different (incorrect) answer.
    # -------------------------------------------------------------------------
    print("[INFO]", product_label, "Computing pixel-level TINV for HDD/CDD")
    hdd_tinv_pixel = ds_daily["HDD20"].groupby("time.year").sum().mean("year")
    cdd_tinv_pixel = ds_daily["CDD20"].groupby("time.year").sum().mean("year")

    # -------------------------------------------------------------------------
    # Main year loop: build all variables (temperature, precip, cross terms)
    # and aggregate to country level in a single batched xa.aggregate call per year.
    # Batching all variables together avoids repeated xagg overhead, and a single
    # .load() computes all numerators at once rather than variable by variable.
    # -------------------------------------------------------------------------
    all_year_rows: list[pd.DataFrame] = []

    for yr in years:
        yr_int = int(yr)
        print(f"[INFO] {product_label} Aggregating {yr_int}")
        ds_y = ds_daily.sel(time=ds_daily.time.dt.year == yr_int)

        ds_sum = xr.Dataset()

        # --- Temperature-derived variables: sum over days in year ---
        for name, da in ds_y.data_vars.items():
            ds_sum[name] = ensure_2d_lat_lon(da, time_agg="sum")

        # --- Precipitation: convert flux to mm/day and sum over days ---
        # PR2 is the sum of squared daily values, capturing within-year
        # intensity/variability consistent with polynomial precip controls
        # in Carleton et al.
        if ds_p is not None and PRDAY_VAR_NAME in ds_p:
            pr_day = ds_p[PRDAY_VAR_NAME].sel(
                time=slice(f"{yr_int}-01-01", f"{yr_int}-12-31")
            )
            if pr_day.sizes.get("time", 0) == 0:
                print(f"[WARN] {product_label} No precip data for year {yr_int}")
            else:
                units = (pr_day.attrs.get("units", "") or "").lower()
                pr_mm_day = (
                    pr_day * 86400.0
                    if "kg m-2 s-1" in units or "kg m**-2 s**-1" in units
                    else pr_day
                )
                ds_sum["PR1"] = pr_mm_day.sum("time", skipna=True)         # mm/year
                ds_sum["PR2"] = (pr_mm_day ** 2).sum("time", skipna=True)  # sum of squared daily precip

        # --- Cross terms: pixel-wise product of annual poly and TINV ---
        # Formed at pixel level here before aggregation so that the
        # population-weighted mean is of the pixel-level product, not
        # the product of two separately aggregated means.
        for i in range(1, 5):
            abv_i = ensure_2d_lat_lon(ds_y[f"ABV{i}"], time_agg="sum")
            blw_i = ensure_2d_lat_lon(ds_y[f"BLW{i}"], time_agg="sum")
            ds_sum[f"ABV{i}_x_CDD"] = abv_i * cdd_tinv_pixel
            ds_sum[f"BLW{i}_x_HDD"] = blw_i * hdd_tinv_pixel

        # --- Batch population-weight all variables and aggregate in one call ---
        # Multiply each gridded field by pixel population (numerator), then load
        # the full dataset at once to avoid repeated per-variable dask compute calls.
        ds_weighted = xr.Dataset({
            name: gpw_on_grid * da.astype("float32").broadcast_like(gpw_on_grid)
            for name, da in ds_sum.data_vars.items()
        }).load()

        with xa.set_options(silent=True):
            agg_all = xa.aggregate(ds_weighted, wm)

        df_year = agg_all.to_dataframe().reset_index()
        idx_col = "shape_id" if "shape_id" in df_year.columns else df_year.columns[0]
        df_year["iso"] = df_year[idx_col].map(shape_to_iso)

        # Divide aggregated numerator by country population to get weighted mean
        df_year = df_year.merge(pop_by_iso.rename("pop_country"), on="iso", how="left")
        for name in ds_sum.data_vars:
            df_year[name] = df_year[name] / df_year["pop_country"]

        df_year["year"] = yr_int
        all_year_rows.append(df_year)

    df = pd.concat(all_year_rows, ignore_index=True)

    # -------------------------------------------------------------------------
    # Rename internal variable names to CIL output convention.
    # Cross term names are included here so everything is renamed in one pass.
    # -------------------------------------------------------------------------
    rename_map = {
        "T1": "temp1_{s}", "T2": "temp2_{s}", "T3": "temp3_{s}", "T4": "temp4_{s}",
        "HDD20": "hdd20_{s}", "CDD20": "cdd20_{s}",
        "ABV1": "polyAbove1_{s}", "ABV2": "polyAbove2_{s}",
        "ABV3": "polyAbove3_{s}", "ABV4": "polyAbove4_{s}",
        "BLW1": "polyBelow1_{s}", "BLW2": "polyBelow2_{s}",
        "BLW3": "polyBelow3_{s}", "BLW4": "polyBelow4_{s}",
        "PR1": "precip1_{s}", "PR2": "precip2_{s}",
        "ABV1_x_CDD": "polyAbove1_x_cdd_{s}", "ABV2_x_CDD": "polyAbove2_x_cdd_{s}",
        "ABV3_x_CDD": "polyAbove3_x_cdd_{s}", "ABV4_x_CDD": "polyAbove4_x_cdd_{s}",
        "BLW1_x_HDD": "polyBelow1_x_hdd_{s}", "BLW2_x_HDD": "polyBelow2_x_hdd_{s}",
        "BLW3_x_HDD": "polyBelow3_x_hdd_{s}", "BLW4_x_HDD": "polyBelow4_x_hdd_{s}",
    }
    rename_map = {k: v.format(s=suffix) for k, v in rename_map.items()}

    keep_cols = ["iso", "year"] + [k for k in rename_map if k in df.columns]
    df = df[keep_cols].rename(columns=rename_map)
    df = df.sort_values(["iso", "year"]).reset_index(drop=True)

    # -------------------------------------------------------------------------
    # Country-level TINV: long-run mean of HDD/CDD averaged across years,
    # derived from the already-aggregated country-year values.
    # These are the country-level heterogeneity moderators used directly
    # in the regression (distinct from the pixel-level TINV used for cross terms above).
    # -------------------------------------------------------------------------
    hdd_col = f"hdd20_{suffix}"
    cdd_col = f"cdd20_{suffix}"

    if hdd_col in df.columns:
        hdd_tinv_country = df.groupby("iso")[hdd_col].mean()
        df[f"hdd20_TINV_{suffix}"] = df["iso"].map(hdd_tinv_country)

    if cdd_col in df.columns:
        cdd_tinv_country = df.groupby("iso")[cdd_col].mean()
        df[f"cdd20_TINV_{suffix}"] = df["iso"].map(cdd_tinv_country)

    print(f"[INFO] {product_label} Writing {out_csv}")
    df.to_csv(out_csv, index=False)
    print(f"[DONE] {product_label}")


# ---------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------


def main() -> None:
    # Loops over all products in car_paths.csv. To restrict to a single product,
    # add a conditional here, e.g.: if product_label != "MERRA2": continue
    car_paths_path = (
        "/user/ab5405/summeraliaclimate/code/energy_uncertainty/car_paths.csv"
    )
    print("[INFO] Reading", car_paths_path)
    car_paths = pd.read_csv(car_paths_path)
    car_paths = car_paths.dropna(subset=["product"])

    for _, row in car_paths.iterrows():
        product_label = str(row["product"]).strip()

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