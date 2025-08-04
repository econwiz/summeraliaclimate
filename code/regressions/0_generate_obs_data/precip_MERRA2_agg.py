#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8
"""
Generate ADM2‐level temperature (MERRA2) and precipitation (GMFD) aggregates for regression,
using the same “one‐loop” logic as your JRA script but with separate weightmaps to avoid grid
mismatch.

- Time‐chunking remains {"time": 1}.
- We compute two weightmaps (tas vs pr) from their native grids.
- In each year, we build one ds_poly containing T1–T4 and PRCP_sum, then aggregate the
  temperature and precipitation pieces separately using their own weightmaps.
- Finally we merge the two ADM2 tables and append the ADM1 long‐run mean temp.
"""

import xagg as xa
import geopandas as gpd
import xarray as xr
import pandas as pd
from pathlib import Path

# ── Config ─────────────────────────────────────────────────────────────────────
product    = "MERRA2"   # label for temperature vars
pr_product = "GMFD"     # label for precipitation vars
BASE       = Path(__file__).resolve().parent
SHAPE      = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
OUTDIR     = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs")
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_CSV    = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"
YEARS      = range(1981, 2014)

# ── Helpers ────────────────────────────────────────────────────────────────────
def open_climate(path: str) -> xr.Dataset:
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path)

def get_precip_var(ds: xr.Dataset) -> str:
    for v in ("prcp", "pr", "precip", "tp"):
        if v in ds.data_vars:
            return v
    raise KeyError(f"No precip var in {ds}")

# ── Paths & shapefile ──────────────────────────────────────────────────────────
df_paths   = pd.read_csv(BASE/"car_paths.csv", dtype=str)
df_paths["product"] = df_paths["product"].str.strip()
row_t      = df_paths.query("product == @product").iloc[0]
row_p      = df_paths.query("product == @pr_product").iloc[0]
tas_path   = row_t["tas_filepath"]
prc_path   = row_p["precip_filepath"]

gdf_adm2   = gpd.read_file(SHAPE, engine="fiona").to_crs("EPSG:4326")
gdf_adm2["geometry"] = gdf_adm2.buffer(0)
gdf_adm2 = gdf_adm2[gdf_adm2.is_valid & ~gdf_adm2.is_empty]
gdf_adm1 = gdf_adm2.dissolve(by=["iso","adm1_id"], as_index=False)

# ── Load & chunk climate data ──────────────────────────────────────────────────
ds_tas = open_climate(tas_path).chunk({"time": 365})
ds_tas = ds_tas.assign(tas=ds_tas.tas - 273.15)
ds_pr  = open_climate(prc_path).chunk({"time": 365})
precip_var = get_precip_var(ds_pr)

# ── Build separate weightmaps ──────────────────────────────────────────────────
sample_t = ds_tas.isel(time=0, drop=True)
sample_p = ds_pr.isel(time=0, drop=True)
wm2_t    = xa.pixel_overlaps(sample_t, gdf_adm2)
wm1_t    = xa.pixel_overlaps(sample_t, gdf_adm1)
wm2_p    = xa.pixel_overlaps(sample_p, gdf_adm2)

# ── ADM1 long‐run mean temp (1981–2010) ────────────────────────────────────────
tmean_lr = ds_tas.sel(time=slice("1981-01-01","2010-12-31")).tas.mean("time")
with xa.set_options(impl="numba", silent=True):
    agg_lr = xa.aggregate(tmean_lr.chunk({"lat":-1,"lon":-1}), wm1_t)
df_lr = (
    agg_lr.to_dataframe().reset_index()
    .rename(columns={"tas": f"lr_tavg_{product}_adm1_avg"})
)[["iso","adm1_id",f"lr_tavg_{product}_adm1_avg"]]

# ── Yearly loop (one big ds_poly, two aggregates) ──────────────────────────────
rows = []
for yr in YEARS:
    print(f"Year {yr}")
    ds_yr = xr.merge([
        ds_tas.sel(time=slice(f"{yr}-01-01",f"{yr}-12-31")),
        ds_pr.sel(time=slice(f"{yr}-01-01",f"{yr}-12-31"))[[precip_var]]
    ]).chunk({"time":1})
    ds_yr = ds_yr.assign(tas=ds_yr.tas - 273.15)

    # build polynomials + precip sum
    ds_poly = xr.merge([
        ds_yr.tas.sum("time").rename("T1_sum"),
        (ds_yr.tas**2).sum("time").rename("T2_sum"),
        (ds_yr.tas**3).sum("time").rename("T3_sum"),
        (ds_yr.tas**4).sum("time").rename("T4_sum"),
        ds_yr[precip_var].sum("time").rename("PRCP_sum"),
    ])

    # aggregate temperature polys
    with xa.set_options(impl="numba", silent=True):
        agg_t = xa.aggregate(ds_poly[["T1_sum","T2_sum","T3_sum","T4_sum"]]
                             .chunk({"lat":-1,"lon":-1}),
                             wm2_t)
    df_t = (
        agg_t.to_dataframe().reset_index()
        .rename(columns={
            "T1_sum":    f"tavg_poly_1_{product}",
            "T2_sum":    f"tavg_poly_2_{product}",
            "T3_sum":    f"tavg_poly_3_{product}",
            "T4_sum":    f"tavg_poly_4_{product}"})
    )

    # aggregate precipitation
    with xa.set_options(impl="numba", silent=True):
        agg_p = xa.aggregate(ds_poly[["PRCP_sum"]].chunk({"lat":-1,"lon":-1}),
                             wm2_p)
    df_p = (
        agg_p.to_dataframe().reset_index()
        .rename(columns={"PRCP_sum": f"prcp_annual_{pr_product}"})
    )
    for p in (1,2,3,4):
        df_p[f"prcp_poly_{p}_{pr_product}"] = df_p[f"prcp_annual_{pr_product}"]**p

    # merge temp & pr, add year & join ADM1 LR temp
    df = (
        pd.merge(df_t, df_p,
                 on=["iso","adm1_id","adm2_id"], how="inner")
          .assign(year=yr, product=product)
          .merge(df_lr, on=["iso","adm1_id"], how="left")
    )
    rows.append(df)

# ── Final combine & write ──────────────────────────────────────────────────────
out = pd.concat(rows, ignore_index=True)
out.to_csv(OUT_CSV, index=False)
print(f"✅ Wrote {OUT_CSV.name} with {len(out)} rows")
