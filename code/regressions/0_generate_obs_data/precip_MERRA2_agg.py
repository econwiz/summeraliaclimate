#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8
"""
This script generates observational temperature and precipitation data at the ADM2 level
from raw climate model files using the xagg package for spatial aggregation.

Steps:
1. Load raw temperature and precipitation data.
2. Preprocess ADM2-level shapefiles and dissolve them to ADM1 level for computing the lr temperature average.
3. Compute long-run average temperature over the entire sample (1981-2013) at the ADM1 level.
4. Loop over each year and compute:
   - Temperature polynomials (sum of daily tas^p for p=1..4)
   - Annual mean precipitation and its polynomials (prcp^p for p=1..4)
5. Spatially aggregate variables to ADM2 using xagg weight maps.
6. Merge in ADM1-level long-run averages for each year.
7. Write results out as a CSV for merging with the mortality panel later.

Inputs:
- car_paths.csv: Contains paths to temperature and precip files for each product.
- mortality_insample_world.shp: Used for spatial aggregation weights.

Output:
- <product>_by_region_year.csv: Contains ADM2-level temperature and precipitation
  variables for each year, formatted for regression prep.

"""
import xagg as xa
import geopandas as gpd
import xarray as xr
import pandas as pd
from pathlib import Path

#select the correct product 
product = "MERRA2"  # Hybrid product name (for output and labels)
BASE    = Path(__file__).resolve().parent

# load the filepaths from the car_paths filedf_paths = pd.read_csv(BASE / "car_paths.csv", dtype=str)
df_paths = pd.read_csv(BASE / "car_paths.csv", dtype=str)
df_paths["product"] = df_paths["product"].str.strip()

 
row_merra = df_paths.loc[df_paths["product"] == "MERRA2"].iloc[0]  # Temperature
row_gmfd  = df_paths.loc[df_paths["product"] == "GMFD"].iloc[0]   # Precipitation

#Read the shape_path from the original Carleton data. This shapefile requires to to_crs command to format it 
#correct as a projection.
shape_path = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

gdf_adm2 = (gpd.read_file(shape_path, engine="fiona").to_crs("EPSG:4326"))
#fix broken geometries, and drop unusable ones.
gdf_adm2["geometry"] = gdf_adm2.buffer(0)
gdf_adm2 = gdf_adm2[gdf_adm2.is_valid & ~gdf_adm2.is_empty]

#Aggregate ADM2 geometries to the ADM1 level by grouping all adm2 polygons that share the same iso and adm1_id into the same 
#polygon. This is for calculating the lr average temp at the adm1 level. Store iso and adm1 as columns rather than the index.
gdf_adm1 = gdf_adm2.dissolve(by=["iso", "adm1_id"], as_index=False)

# Helper function to open the climate data, which may be a netCDF or a zarr.
def open_climate(path):
    if path.endswith(".zarr"):
        return xr.open_zarr(path, consolidated=False)
    return xr.open_dataset(path)

#Open the climate and precipitation files, and applies chunking on the time dimension, so that only one time step is read into memory at a time.
ds_tas = open_climate(row_merra["tas_filepath"]).chunk({"time": 1})
ds_pr  = open_climate(row_gmfd["precip_filepath"]).chunk({"time": 1})


#find what the precip variable is called
precip_var = None
for cand in ("prcp", "pr", "precip"):
    if cand in ds_pr.data_vars:
        precip_var = cand
        break
if precip_var is None:
    raise KeyError(f"no recognized precip var in {row_gmfd['precip_filepath']}")

#merge the temperature data and the precipitation data, again chunking along the time dimension. Convert the temps from K to C.
ds_full = xr.merge([ds_tas, ds_pr[[precip_var]]]).chunk({"time": 1})
ds_full = ds_full.assign(tas = ds_full.tas - 273.15)

#Build the weightmaps. To reduce memory usage, it grabs a single timeslice (because the spatial boundaries are constant over time), therefore dissolving it to a latxlon grid. Then, calculate two weightmaps, one at the adm2 level from the original shapefile and one with the dissolved shapefile.
sample   = ds_full.isel(time=0, drop=True)
wm_adm2  = xa.pixel_overlaps(sample, gdf_adm2)
wm_adm1  = xa.pixel_overlaps(sample, gdf_adm1)

# Compute the lr average temperature for each adm1 region over the entire sample period, to use as a control later. 
clim_start, clim_end = 1981, 2010
ds_clim  = ds_full.sel(time=slice(f"{clim_start}-01-01", f"{clim_end}-12-31"))
tmean_lr = ds_clim.tas.mean("time").chunk({"lat": -1, "lon": -1}).rename("Tmean")

#aggregate the temperature field for the lr average temp on the adm1 level, produces an xarray dataset with one row per adm1 region.
with xa.set_options(impl="numba", silent=True):
    agg_lr = xa.aggregate(tmean_lr, wm_adm1)

#Converts the aggregated temperature dataset to a pandas df, select only relevant columns.
df_lr = (
    agg_lr.to_dataframe().reset_index()[['iso','adm1_id','Tmean']]
        .rename(columns={'Tmean': f"lr_tavg_{product.upper()}_adm1_avg"})
)

#Set the name of the output and select years.
years   = range(1981, 2014)
OUTDIR = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs")
out_csv = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"

#for each year in the sample list, calculate the temperature polynomials of the sum of daily temp daily over the year. First, get the year by selecting jan 1 - dec 31 and sum the daily temps from that period for each year. Also sum daily precip data over the year.
rows = []
for yr in years:
    print(f"  Year {yr}")
    ds_yr = ds_full.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))

    # Temperature polynomials
    t1 = ds_yr.tas.sum("time").rename("T1_sum")
    t2 = (ds_yr.tas ** 2).sum("time").rename("T2_sum")
    t3 = (ds_yr.tas ** 3).sum("time").rename("T3_sum")
    t4 = (ds_yr.tas ** 4).sum("time").rename("T4_sum")

    # Precipitation sum
    pr_sum = ds_yr[precip_var].sum("time").rename("PRCP_sum")

    ds_poly = xr.merge([t1, t2, t3, t4, pr_sum]).chunk({"lat": -1, "lon": -1})
    
    #aggregate the polynomial data as the adm2 level, i.e. to the original shapefile level.
    with xa.set_options(impl="numba", silent=True):
        agg = xa.aggregate(ds_poly, wm_adm2)
        
    #transform the data to a pandas df.
    df_reg = (
        agg.to_dataframe().reset_index()
            .rename(columns={
                "T1_sum":    f"tavg_poly_1_{product.upper()}",
                "T2_sum":    f"tavg_poly_2_{product.upper()}",
                "T3_sum":    f"tavg_poly_3_{product.upper()}",
                "T4_sum":    f"tavg_poly_4_{product.upper()}",
                "PRCP_sum": f"prcp_annual_{product.upper()}"
            })
            .assign(year=yr, product=product)
    )

    # Precipitation polynomials
    for p in [1, 2, 3, 4]:
        df_reg[f"prcp_poly_{p}_{product.upper()}"] = df_reg[f"prcp_annual_{product.upper()}"] ** p

    #merge the precip data and the temp data for each year on iso and adm1 id, which they share. 
    df_reg = df_reg.merge(df_lr, on=['iso','adm1_id'], how='left')
    rows.append(df_reg)

# Combine yearly dataframes into one master dataframe, replacing any missing values with 0. Write the full output to CSV for use in downstream regression scripts.
prod_df = pd.concat(rows, ignore_index=True)
print("Writing output to:", out_csv)
prod_df.to_csv(out_csv, index=False)
print(f"↳ Wrote {out_csv.name} ({len(prod_df)} rows)")
print("✅ Done!")
