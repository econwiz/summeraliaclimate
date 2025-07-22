#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8

import xagg as xa
import geopandas as gpd 
import xarray as xr
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent

df_paths = pd.read_csv(BASE / "car_paths.csv", dtype=str)

shape_path = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

gdf = gpd.read_file(shape_path, engine="fiona")
gdf["geometry"] = gdf.buffer(0)
gdf = gdf[gdf.is_valid & (~gdf.is_empty)]

def open_climate(path):
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    else:
        return xr.open_dataset(path)
    
years = range(1981,2014)
OUTDIR = Path("/user/ab5405/summeraliaclimate/code/obs_csvs")
OUTDIR.mkdir(parents=True, exist_ok=True)

for _, row in df_paths.iterrows():
    product = row["product"].strip()
    
    out_csv = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"

    # skip if it already exists
    if out_csv.exists():
        print(f"Skipping {product}: {out_csv.name} already exists")
        continue

    path    = row["filepath"]
    print(f"Testing product: {product}")

    ds_full = open_climate(path).chunk({"time": 1})
    ds_full = ds_full.assign(tas=ds_full.tas - 273.15)
    
    sample = ds_full.isel(time=0, drop=True)
    wm     = xa.pixel_overlaps(sample, gdf)
    rows_for_prod = []

    for yr in years:
        print(f"  Year {yr}")
        
        ds_yr = ds_full.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))
        
        t1    = ds_yr.tas.sum("time").rename("T1_sum")
        t2    = (ds_yr.tas ** 2).sum("time").rename("T2_sum")
        t3    = (ds_yr.tas ** 3).sum("time").rename("T3_sum")
        t4    = (ds_yr.tas ** 4).sum("time").rename("T4_sum")
        tmean = ds_yr.tas.mean("time").rename("Tmean")

        ds_poly = xr.merge([t1, t2, t3, t4, tmean])
        ds_poly = ds_poly.chunk({"lat": -1, "lon": -1})
        
        with xa.set_options(impl="numba", silent=True):
            agg = xa.aggregate(ds_poly,wm)
            
        df_reg = agg.to_dataframe().reset_index()
        safe = product.upper().replace("-", "_")  
        
        df_reg = (agg
          .to_dataframe()
          .reset_index()
          .rename(columns={
             "T1_sum": f"tavg_poly_1_{product.upper()}",
             "T2_sum": f"tavg_poly_2_{product.upper()}",
             "T3_sum": f"tavg_poly_3_{product.upper()}",
             "T4_sum": f"tavg_poly_4_{product.upper()}",
             "Tmean" : f"lr_tavg_{product.upper()}_adm1_avg",
          })
          .assign(product=product, year=yr)
        )
        rows_for_prod.append(df_reg)
        
    prod_df = pd.concat(rows_for_prod, ignore_index=True)
    prod_df = prod_df.fillna(0)
    prod_df.to_csv(out_csv, index=False)
    print(f"   ↳ Wrote {out_csv.name} ({len(prod_df)} rows)")

print("✅ All products complete!")
