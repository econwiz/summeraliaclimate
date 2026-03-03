#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8

import xagg as xa
import geopandas as gpd
import xarray as xr
import pandas as pd
from pathlib import Path

product = "GMFD"   # Change to ERA5-025, GMFD, etc. as needed
BASE    = Path(__file__).resolve().parent

# ─── 1) Load climate file paths ────────────────────────────────────────────────
df_paths = pd.read_csv(BASE / "car_paths.csv", dtype=str)
df_paths["product"] = df_paths["product"].str.strip()
row = df_paths.loc[df_paths["product"] == product].iloc[0]

# ─── 2) Read ADM2 shapefile and derive ADM1 geometry ───────────────────────────
shape_path = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
gdf_adm2   = gpd.read_file(shape_path, engine="fiona").to_crs("EPSG:4326")
gdf_adm2["geometry"] = gdf_adm2.buffer(0)
gdf_adm2   = gdf_adm2[gdf_adm2.is_valid & ~gdf_adm2.is_empty]
gdf_adm1   = gdf_adm2.dissolve(by=["iso", "adm1_id"], as_index=False)

# ─── 3) Helper: open climate data ─────────────────────────────────────────────
def open_climate(path):
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path)

print(f"Processing {product}")

# ─── 4) Load temperature data and convert to °C ───────────────────────────────
ds_tas  = open_climate(row["tas_filepath"])
ds_full = ds_tas.chunk({"time": 1})
ds_full = ds_full.assign(tas=ds_full.tas - 273.15)

# ─── 5) Build weight maps for ADM2 & ADM1 ─────────────────────────────────────
sample   = ds_full.isel(time=0, drop=True)
wm_adm2  = xa.pixel_overlaps(sample, gdf_adm2)
wm_adm1  = xa.pixel_overlaps(sample, gdf_adm1)

# ─── 5b) Compute long-run ADM1 mean (once) ────────────────────────────────────
clim_start, clim_end = 1981, 2010  # adjust to match your climate window
ds_clim  = ds_full.sel(time=slice(f"{clim_start}-01-01", f"{clim_end}-12-31"))
tmean_lr = ds_clim.tas.mean("time").chunk({"lat": -1, "lon": -1}).rename("Tmean")

with xa.set_options(impl="numba", silent=True):
    agg1_lr = xa.aggregate(tmean_lr, wm_adm1)

df1_lr = (
    agg1_lr.to_dataframe().reset_index()[['iso', 'adm1_id', 'Tmean']]
          .rename(columns={'Tmean': f"lr_tavg_{product.upper()}_adm1_avg"})
)

# ─── 6) Loop over years and aggregate ─────────────────────────────────────────
years   = range(1981, 2014)
OUTDIR = BASE / "generate_obs_data" / "obs_csvs"      # <─ NEW path
OUTDIR.mkdir(parents=True, exist_ok=True)
out_csv = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"
rows    = []

for yr in years:
    print(f"  Year {yr}")
    ds_yr = ds_full.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))

    # temperature polynomial sums
    t1    = ds_yr.tas.sum("time").rename("T1_sum")
    t2    = (ds_yr.tas**2).sum("time").rename("T2_sum")
    t3    = (ds_yr.tas**3).sum("time").rename("T3_sum")
    t4    = (ds_yr.tas**4).sum("time").rename("T4_sum")
    tmean = ds_yr.tas.mean("time").rename("Tmean")

    # merge into one Dataset
    ds_poly = xr.merge([t1, t2, t3, t4, tmean]).chunk({"lat": -1, "lon": -1})

    with xa.set_options(impl="numba", silent=True):
        # aggregate at ADM2
        agg2 = xa.aggregate(ds_poly, wm_adm2)

    # to pandas: ADM2-level
    df2 = (
        agg2.to_dataframe().reset_index()
            .rename(columns={
                "T1_sum": f"tavg_poly_1_{product.upper()}",
                "T2_sum": f"tavg_poly_2_{product.upper()}",
                "T3_sum": f"tavg_poly_3_{product.upper()}",
                "T4_sum": f"tavg_poly_4_{product.upper()}",
                "Tmean":  f"tmean_{product.upper()}_adm2_year"
            })
            .assign(year=yr, product=product)
    )

    # merge the long-run ADM1 mean (already computed)
    df_reg = df2.merge(df1_lr, on=['iso', 'adm1_id'], how='left')
    rows.append(df_reg)

# ─── 7) Finalize and write ────────────────────────────────────────────────────
prod_df = pd.concat(rows, ignore_index=True).fillna(0)
prod_df.to_csv(out_csv, index=False)
print(f"↳ Wrote {out_csv.name} ({len(prod_df)} rows)")
print("✅ Done!")
