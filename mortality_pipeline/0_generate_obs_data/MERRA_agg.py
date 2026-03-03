#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8

"""
This script generates observational temperature data at the ADM2 level
from raw climate model files using the xagg package for spatial aggregation.

Steps:
1. Load raw temperature data.
2. Preprocess ADM2-level shapefiles and dissolve them to ADM1 level for computing the lr temperature average.
3. Compute long-run average temperature over the entire sample (1981-2013) at the ADM1 level.
4. Loop over each year and compute:
   - Temperature polynomials (sum of daily tas^p for p=1..4)
   # - Annual precipitation (commented out)
5. Spatially aggregate variables to ADM2 using xagg weight maps.
6. Merge in ADM1-level long-run averages for each year.
7. Write results out as a CSV for merging with the mortality panel later.

Inputs:
- car_paths.csv: Contains paths to temperature (and precip) files for each product.
- mortality_insample_world.shp: Used for spatial aggregation weights.

Output:
- <product>_by_region_year.csv: Contains ADM2-level temperature variables for each year, formatted for regression prep.
"""

import xagg as xa
import geopandas as gpd
import xarray as xr
import pandas as pd
from pathlib import Path

#select the correct product 
product = "MERRA2"
BASE    = Path(__file__).resolve().parent

# load the filepaths from the car_paths file
df_paths = pd.read_csv(BASE / "car_paths.csv", dtype=str)
df_paths["product"] = df_paths["product"].str.strip()
row        = df_paths.loc[df_paths["product"] == product].iloc[0]

#Read the shape_path from the original Carleton data. This shapefile requires to to_crs command to format it 
#correct as a projection. 
shape_path = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"
gdf_adm2  = gpd.read_file(shape_path, engine="fiona").to_crs("EPSG:4326")
#fix broken geometries, and drop unusable ones.
gdf_adm2["geometry"] = gdf_adm2.buffer(0)
gdf_adm2  = gdf_adm2[gdf_adm2.is_valid & ~gdf_adm2.is_empty]
#Aggregate ADM2 geometries to the ADM1 level by grouping all adm2 polygons that share the same iso and adm1_id into the same 
#polygon. This is for calculating the lr average temp at the adm1 level. Store iso and adm1 as columns rather than the index.
gdf_adm1  = gdf_adm2.dissolve(by=["iso","adm1_id"], as_index=False)

# Helper function to open the climate data, which may be a netCDF or a zarr.
def open_climate(path):
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path)

print(f"Processing {product}")

# Open climate files; chunk time by 1 day
ds_tas = open_climate(row["tas_filepath"]).chunk({"time": 365})
# ds_pr  = open_climate(row["precip_filepath"]).chunk({"time": 1}) 
# precip_var = next((v for v in ("prcp","pr","precip","tp") if v in ds_pr.data_vars), None)
# if precip_var is None:
#     raise KeyError(f"No precip var in {row['precip_filepath']}; got {list(ds_pr.data_vars)}")

# Merge (temperature only) and convert K -> C
# ds_full = xr.merge([ds_tas, ds_pr[[precip_var]]]).chunk({"time": 1})  # <-- commented out
ds_full = ds_tas.chunk({"time": 1})
ds_full = ds_full.assign(tas = ds_full.tas - 273.15)

# Build weightmaps (from a single timeslice)
sample   = ds_full.isel(time=0, drop=True)
wm_adm2  = xa.pixel_overlaps(sample, gdf_adm2)
wm_adm1  = xa.pixel_overlaps(sample, gdf_adm1)

# Compute lr average temperature for ADM1 (1981–2010)
clim_start, clim_end = 1981, 2010
ds_clim  = ds_full.sel(time=slice(f"{clim_start}-01-01", f"{clim_end}-12-31"))
tmean_lr = ds_clim.tas.mean("time").chunk({"lat": -1, "lon": -1}).rename("Tmean")

# Aggregate lr temp to ADM1
with xa.set_options(impl="numba", silent=True):
    agg_lr = xa.aggregate(tmean_lr, wm_adm1)

# Convert ADM1 lr temp to DataFrame
df_lr = (
    agg_lr
    .to_dataframe()
    .reset_index()[['iso','adm1_id','Tmean']]
    .rename(columns={'Tmean': f"lr_tavg_{product.upper()}_adm1_avg"})
)

# Years & output path
years   = range(1981, 2014)
OUTDIR = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs")
out_csv = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"

# Yearly aggregation (temperature only; precip commented out)
rows = []
for yr in years:
    print(f"  Year {yr}")
    ds_yr = ds_full.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))

    # Temperature polynomials (sum over days)
    t1 = ds_yr.tas.sum("time").rename("T1_sum")
    t2 = (ds_yr.tas**2).sum("time").rename("T2_sum")
    t3 = (ds_yr.tas**3).sum("time").rename("T3_sum")
    t4 = (ds_yr.tas**4).sum("time").rename("T4_sum")

    # Precipitation (commented out)
    # pr_sum = ds_yr[precip_var].sum("time").rename("PRCP_sum")

    # Merge fields and single-chunk spatial dims
    # ds_poly = xr.merge([t1, t2, t3, t4, pr_sum]).chunk({"lat": -1, "lon": -1})  # <-- commented precip
    ds_poly = xr.merge([t1, t2, t3, t4]).chunk({"lat": -1, "lon": -1})

    # Aggregate to ADM2
    with xa.set_options(impl="numba", silent=True):
        agg2 = xa.aggregate(ds_poly, wm_adm2)

    # To DataFrame & rename columns
    df = (
        agg2
        .to_dataframe()
        .reset_index()
        .rename(columns={
            "T1_sum":    f"tavg_poly_1_{product.upper()}",
            "T2_sum":    f"tavg_poly_2_{product.upper()}",
            "T3_sum":    f"tavg_poly_3_{product.upper()}",
            "T4_sum":    f"tavg_poly_4_{product.upper()}",
            # "PRCP_sum": f"prcp_annual_{product.upper()}",   # <-- commented out
        })
        .assign(year=yr, product=product)
    )

    # Precip polynomials (commented out)
    # for p in (1,2,3,4):
    #     df[f"prcp_poly_{p}_{product.upper()}"] = df[f"prcp_annual_{product.upper()}"] ** p

    # Merge ADM1 lr temp
    df = df.merge(df_lr, on=['iso','adm1_id'], how='left')
    rows.append(df)

# Combine years & write CSV
out = pd.concat(rows, ignore_index=True)
out.to_csv(out_csv, index=False)
print(f"↳ Wrote {out_csv.name} with {len(out)} rows")
