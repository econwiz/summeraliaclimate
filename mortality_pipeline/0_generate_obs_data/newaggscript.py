#!/usr/bin/env python
# coding: utf-8
"""
Generate ADM2-by-year climate inputs for ALL products in car_paths.csv.

Matches Kevin's approach:
- Geometry: make_valid() only (no buffer(0)); WGS84; drop empties.
- Years: align to the mortality panel years.
- Temperature (°C): yearly sums of tas^k, k=1..4, at ADM2.
- Precipitation: yearly sums of pr^p, p=1..4, at ADM2 (Kevin-style).
- Long-run mean T (1981–2010) at ADM1.
- float32 everywhere to keep files light.
- Chunking tuned for performance and stable graphs.
- If a product has no precip (e.g., MERRA2), keep precip columns with zeros.
"""

import warnings
warnings.filterwarnings("ignore")

import time
from pathlib import Path
import pandas as pd
import geopandas as gpd
import xarray as xr
import xagg as xa
from shapely.validation import make_valid

BASE = Path(__file__).resolve().parent
CAR_PATHS = BASE / "car_paths.csv"

# Original Carleton ADM2 shapefile
SHAPE_PATH = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

# Mortality panel (to align years)
PANEL_DTA = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"

# Output directory
OUTDIR = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs")
OUTDIR.mkdir(parents=True, exist_ok=True)

# Long-run climatology window
CLIM_START, CLIM_END = 1981, 2010

# -------------------- Helpers --------------------
def open_climate(path: str) -> xr.Dataset:
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path)


def repair_geoms_like_kevin(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf = gdf.copy()
    gdf["geometry"] = gdf["geometry"].make_valid()
    gdf = gdf.to_crs("EPSG:4326")


def detect_precip_var(ds: xr.Dataset) -> str | None:
    return next((v for v in ("pr", "prcp", "precip", "tp") if v in ds.data_vars), None)


def normalize_precip_units(da: xr.DataArray, name: str) -> xr.DataArray:
    """Normalize precip to mm/day so daily sums are meaningful (Kevin’s CSVs)."""
    units = str((da.attrs or {}).get("units", "")).lower()
    out = da.astype("float32")
    if name == "pr" and (("kg" in units and "s" in units) or units in ("kg m**-2 s**-1", "kg/m^2/s", "kg m-2 s-1")):
        out = out * 86400.0
        out.attrs["units"] = "mm/day"
    return out


def ensure_float32(ds: xr.Dataset, names: list[str]) -> xr.Dataset:
    for n in names:
        if n in ds:
            ds[n] = ds[n].astype("float32")
    return ds


def get_years_from_panel() -> range:
    """Mirror Kevin’s years by reading Stata panel; fallback to 1981–2013."""
    try:
        import pyreadstat
        df, _ = pyreadstat.read_dta(PANEL_DTA, usecols=["year"])
        years = pd.Series(df["year"]).dropna().astype(int).unique()
        if len(years) > 0:
            return range(int(years.min()), int(years.max()) + 1)
    except Exception:
        pass
    return range(1981, 2014)


# -------------------- Main --------------------
def main():
    YEARS = get_years_from_panel()
    print(f"Years to compute (aligned): {YEARS.start}–{YEARS.stop-1}")

    # Read product paths
    df_paths = pd.read_csv(CAR_PATHS, dtype=str)
    df_paths["product"] = df_paths["product"].str.strip()

    # Load & prepare shapefiles (ADM2 + ADM1)
    gdf_adm2 = repair_geoms_like_kevin(gpd.read_file(SHAPE_PATH, engine="fiona"))
    gdf_adm1 = gdf_adm2.dissolve(by=["iso", "adm1_id"], as_index=False)

    for _, row in df_paths.iterrows():
        t0 = time.time()
        product = row["product"]
        tag = product.replace("-", "_").upper()
        print(f"\n=== {product} ===")

        # Open datasets with stable chunking (avoid time=1 graphs)
        ds_tas = open_climate(row["tas_filepath"]).chunk({"lat": -1, "lon": -1, "time": 30})

        pr_ok = True
        pvar = None
        ds_pr = None
        try:
            ds_pr = open_climate(row["precip_filepath"]).chunk({"lat": -1, "lon": -1, "time": 30})
            pvar = detect_precip_var(ds_pr)
            pr_ok = pvar is not None
        except Exception:
            pr_ok = False

        if pr_ok:
            ds_full = xr.merge([ds_tas[["tas"]], ds_pr[[pvar]]]).chunk({"time": 30})
            ds_full = ds_full.assign(
                tas=(ds_full.tas - 273.15).astype("float32"),
                **{pvar: normalize_precip_units(ds_full[pvar], pvar)}
            )
            ds_full = ensure_float32(ds_full, ["tas", pvar])
        else:
            # temperature-only pipeline; we’ll keep precip columns as zeros later
            ds_full = ds_tas[["tas"]].assign(tas=(ds_tas.tas - 273.15).astype("float32")).chunk({"time": 30})
            ds_full = ensure_float32(ds_full, ["tas"])

        # Build weight maps from a sample slice (Kevin uses xa.pixel_overlaps)
        sample = ds_full.isel(time=0, drop=True)
        with xa.set_options(silent=True):
            wm_adm2 = xa.pixel_overlaps(sample, gdf_adm2)
            wm_adm1 = xa.pixel_overlaps(sample, gdf_adm1)

        # Long-run mean T (1981–2010) at ADM1
        ds_clim = ds_full.sel(time=slice(f"{CLIM_START}-01-01", f"{CLIM_END}-12-31"))
        tmean_lr = ds_clim.tas.mean("time").chunk({"lat": -1, "lon": -1}).rename("Tmean")
        with xa.set_options(impl="numba", silent=True):
            agg_lr = xa.aggregate(tmean_lr, wm_adm1)
        df_lr = (
            agg_lr.to_dataframe()
            .reset_index()[["iso", "adm1_id", "Tmean"]]
            .rename(columns={"Tmean": f"lr_tavg_{tag}_adm1_avg"})
        )

        # Per-year loop (ADM2 polys)
        rows = []
        for yr in YEARS:
            t_year = time.time()
            ds_yr = ds_full.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))

            # Temperature polynomials: sum over days of (tas^k)
            t1 = ds_yr.tas.sum("time").rename("T1_sum")
            t2 = (ds_yr.tas ** 2).sum("time").rename("T2_sum")
            t3 = (ds_yr.tas ** 3).sum("time").rename("T3_sum")
            t4 = (ds_yr.tas ** 4).sum("time").rename("T4_sum")

            if pr_ok:
                pr1 = ds_yr[pvar].sum("time").rename("PR1_sum")
                pr2 = (ds_yr[pvar] ** 2).sum("time").rename("PR2_sum")
                pr3 = (ds_yr[pvar] ** 3).sum("time").rename("PR3_sum")
                pr4 = (ds_yr[pvar] ** 4).sum("time").rename("PR4_sum")
            else:
                # keep schema identical with zeros if precip missing (e.g., MERRA2)
                zeros = xr.zeros_like(t1, dtype="float32")
                pr1 = zeros.rename("PR1_sum")
                pr2 = zeros.rename("PR2_sum")
                pr3 = zeros.rename("PR3_sum")
                pr4 = zeros.rename("PR4_sum")

            ds_poly = xr.merge([t1, t2, t3, t4, pr1, pr2, pr3, pr4]).chunk({"lat": -1, "lon": -1})

            # Aggregate to ADM2
            with xa.set_options(impl="numba", silent=True):
                agg2 = xa.aggregate(ds_poly, wm_adm2)

            # Tidy to DataFrame
            df = (
                agg2.to_dataframe()
                .reset_index()
                .rename(
                    columns={
                        "T1_sum": "tavg_poly_1",
                        "T2_sum": "tavg_poly_2",
                        "T3_sum": "tavg_poly_3",
                        "T4_sum": "tavg_poly_4",
                        "PR1_sum": "prcp_poly_1",
                        "PR2_sum": "prcp_poly_2",
                        "PR3_sum": "prcp_poly_3",
                        "PR4_sum": "prcp_poly_4",
                    }
                )
                .assign(year=yr, product=product)
            )

            # Add product suffix (ERA5_025 etc.) to poly columns
            ren = {
                c: f"{c}_{tag}"
                for c in [
                    "tavg_poly_1",
                    "tavg_poly_2",
                    "tavg_poly_3",
                    "tavg_poly_4",
                    "prcp_poly_1",
                    "prcp_poly_2",
                    "prcp_poly_3",
                    "prcp_poly_4",
                ]
            }
            df = df.rename(columns=ren)

            # Merge LR Tmean at ADM1
            df = df.merge(df_lr, on=["iso", "adm1_id"], how="left")

            rows.append(df)
            print(f"  {yr}: done in {time.time() - t_year:.1f}s")

        # Write output
        out = pd.concat(rows, ignore_index=True)
        out_csv = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"
        out.to_csv(out_csv, index=False)
        print(f"↳ Wrote {out_csv} with {len(out)} rows in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
