#!/apps/anaconda3/bin/python
# coding: utf-8

import xagg as xa
import geopandas as gpd 
import xarray as xr
import pandas as pd
import pyreadstat

df_paths = pd.read_csv("car_paths.csv", dtype=str)

shape_path = "/shared/share_hle/data/aux_data/geo_data/impact-region.shp"
gdf = gpd.read_file(shape_path)

xa.set_options(impl="numba", silent = True)

def open_climate(path):
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    else:
        return xr.open_dataset(path)
    
all_rows = []
years = range(1957,2014)

for _, row in df_paths.iterrows():
    product = row["product"]
    path    = row["filepath"]
    print(f"Processing {product}")

    ds_all = open_climate(path).chunk({"time": 30})
    wm      = xa.pixel_overlaps(ds_grid, gdf)

    for yr in years:
        print(f"  Year {yr}")
        ds_yr = ds_all.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))

        ds_poly = xr.Dataset({
            "T1_sum": (ds_yr.tas).sum("time"),
            "T2_sum": (ds_yr.tas**2).sum("time"),
            "T3_sum": (ds_yr.tas**3).sum("time"),
            "T4_sum": (ds_yr.tas**4).sum("time"),
            "Tmean" :  ds_yr.tas.mean("time"),
        })

        agg    = xa.aggregate(ds_poly, wm)
        df_reg = agg.to_dataframe().reset_index()
        safe = product.upper().replace("-", "_")   # e.g. "ERA5-025" → "ERA5_025"

        df_reg = df_reg.rename(columns={
            "T1_sum": f"tavg_poly_1_{safe}",
            "T2_sum": f"tavg_poly_2_{safe}",
            "T3_sum": f"tavg_poly_3_{safe}",
            "T4_sum": f"tavg_poly_4_{safe}",
            "Tmean" : f"lr_tavg_{safe}_adm1_avg",   # or whatever your prep_data.do expects
        })
        df_reg["product"] = product
        df_reg["year"]    = yr
        all_rows.append(df_reg)

big = pd.concat(all_rows, ignore_index=True)

big = big.rename(columns={"region_id":"region"})
if "time" in big.columns:
    big = big.drop(columns="time")
cols = ["product","region","year"] + [c for c in big.columns if c not in ("product","region","year")]
big = big[cols]

out_path = "/mnt/data/all_products_by_region_year.dta"
big.to_stata(out_path, write_index=False, version=118)
print(f"Wrote {out_path}")
