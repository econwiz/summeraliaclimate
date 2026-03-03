#!/user/ab5405/.conda/envs/hle_iv/bin/python
# coding: utf-8
"""
Generate ADM2-by-year climate inputs aligned to the mortality panel (Kevin parity, simplified).

- EU collapse includes FRA, LIE, TUR for BOTH panel & shapes.
- Build join-only mixed keys (numeric → integer-like strings; non-numeric → trimmed string).
- Keep empties; drop only invalid; dissolve to PANEL ids.
- Hard-gate climate to panel years; tas→°C; precip→mm/day BEFORE powers.
- Long-run T at ADM1: mean of yearly means over 1981–2010 (fallback to panel years).
- FINAL STEP: panel-anchored universe left-join so zero-overlap polygons remain.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import xagg as xa

# ---------------- user-editable constants (kept simple) ----------------
product    = "MERRA2"  # e.g., "GMFD", "JRA_3Q", "ERA5-025", "MERRA2"
BASE       = Path(__file__).resolve().parent

CAR_PATHS  = BASE / "car_paths.csv"  # must map product → tas_filepath, precip_filepath (optional)
PANEL_DTA  = "/shared/share_hle/data/aux_data/global_mortality_panel_public.dta"
SHAPEFILE  = "/shared/share_hle/data/1_estimation/3_regions/insample_shp/mortality_insample_world.shp"

OUTDIR     = Path("/user/ab5405/summeraliaclimate/code/regressions/0_generate_obs_data/obs_csvs")
OUTDIR.mkdir(parents=True, exist_ok=True)
out_csv    = OUTDIR / f"{product.replace('-', '_')}_by_region_year.csv"

# ---------------- helpers ----------------
EU_ISO = [
    "AUT","BEL","BGR","CHE","CYP","CZE","DEU","DNK","ESP",
    "EST","FRA","FIN","GBR","GRC","HRV","HUN","IRL","ISL","ITA",
    "LIE","LTU","LUX","LVA","MKD","MLT","MNE","NLD","NOR",
    "POL","PRT","ROU","SVK","SVN","SWE","TUR"
]

def collapse_eu_iso(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.upper()
    return s.where(~s.isin(EU_ISO), "EU")

def mixed_key(s: pd.Series) -> pd.Series:
    """
    Kevin-style join key:
      - If purely numeric (e.g., '001'), coerce to integer-like string ('1').
      - Otherwise keep the trimmed string (e.g., '01A' stays '01A').
    """
    s_str = s.astype(str).str.strip()
    is_numeric = s_str.str.fullmatch(r"\d+")
    # strip leading zeros only for purely numeric values
    s_num = s_str.where(~is_numeric, s_str.str.lstrip("0").replace({"": "0"}))
    return s_num

def open_any(path):
    p = str(path)
    if p.endswith(".zarr"):
        return xr.open_zarr(p, consolidated=False)
    return xr.open_dataset(p)

def detect_precip_var(ds):
    for v in ("pr", "prcorr", "prcp", "precip", "tp"):
        if v in ds.data_vars:
            return v
    return None

def product_tag(name: str) -> str:
    return name.replace("-", "_").upper()

def normalize_precip_to_mm_per_day(ds: xr.Dataset, var: str, product: str) -> xr.Dataset:
    """
    Convert precip variable to mm/day before polynomial powers.
      - pr/prcorr typically kg m^-2 s^-1 → ×86400
      - tp (ERA5) meters → ×1000  (daily accumulation ⇒ mm/day)
    """
    arr = ds[var]
    units = (arr.attrs.get("units") or "").lower().strip()
    factor = None

    if units in {"kg m**-2 s**-1", "kg m-2 s-1", "kg/m^2/s"}:
        factor = 86400.0
    elif units in {"m", "meter", "metre", "meters"}:
        factor = 1000.0
    elif var in {"pr", "prcorr"}:
        factor = 86400.0
    elif var == "tp":
        factor = 1000.0
    else:
        # fallback by product heuristic
        factor = 1000.0 if product.upper().startswith("ERA5") and var == "tp" else 86400.0

    out = ds.copy()
    out[var] = arr * factor
    out[var].attrs["units"] = "mm/day"
    return out

TAG = product_tag(product)

# ---------------- read inputs ----------------
paths = pd.read_csv(CAR_PATHS, dtype=str)
paths["product"] = paths["product"].str.strip()
if not (paths["product"] == product).any():
    raise ValueError(f"Product '{product}' not found in {CAR_PATHS}")
row = paths.loc[paths["product"] == product].iloc[0]

# Panel (canonical IDs kept; build *_k for joining only)
panel = pd.read_stata(PANEL_DTA)
panel["iso"]  = collapse_eu_iso(panel["iso"])
panel["year"] = pd.to_numeric(panel["year"], errors="coerce").astype("Int64")
panel = panel.dropna(subset=["year"]).copy()
panel["year"] = panel["year"].astype(int)

panel["adm1_k"] = mixed_key(panel["adm1_id"])
panel["adm2_k"] = mixed_key(panel["adm2_id"])

keys = panel[["iso","adm1_id","adm2_id","year","adm1_k","adm2_k"]].drop_duplicates()
spatial_keys = keys[["iso","adm1_id","adm2_id","adm1_k","adm2_k"]].drop_duplicates()
years_insample = np.sort(keys["year"].unique())

print(f"[INFO] Panel years: {years_insample[0]}–{years_insample[-1]}  (n={len(years_insample)})")
print(f"[INFO] Distinct spatial keys in panel: {len(spatial_keys)}")

# Shapefile (keep empties; drop only invalid; dissolve to panel IDs)
gdf_full = gpd.read_file(SHAPEFILE, engine="fiona")
gdf_full["geometry"] = gdf_full["geometry"].make_valid()
gdf_full = gdf_full.to_crs("EPSG:4326")

use_cols = ["iso","adm1_id","adm2_id","geometry"]
missing = [c for c in use_cols if c not in gdf_full.columns]
if missing:
    raise ValueError(f"Shapefile missing required columns: {missing}")

gdf = gdf_full[use_cols].copy()
gdf["iso"] = collapse_eu_iso(gdf["iso"])
gdf = gdf[gdf.is_valid]  # empties are allowed; only invalid are dropped

gdf["adm1_k"] = mixed_key(gdf["adm1_id"])
gdf["adm2_k"] = mixed_key(gdf["adm2_id"])

left  = gdf[["iso","adm1_k","adm2_k","geometry"]].copy()
right = spatial_keys[["iso","adm1_k","adm2_k","adm1_id","adm2_id"]].drop_duplicates()

gdf_sub = left.merge(right, on=["iso","adm1_k","adm2_k"], how="inner", validate="m:1")
gdf_sub = gpd.GeoDataFrame(gdf_sub, geometry="geometry", crs=gdf.crs)
gdf_sub = gdf_sub.dissolve(by=["iso","adm1_id","adm2_id"], as_index=False)
gdf_adm1 = gdf_sub.dissolve(by=["iso","adm1_id"], as_index=False)

print(f"[INFO] Shapes aligned: {len(gdf_sub)} ADM2, {len(gdf_adm1)} ADM1")
print(f"[CHECK] IND ADM2 polygons present: {(gdf_sub['iso']=='IND').sum()}")

# Climate (hard-gate years; tas→°C; precip→mm/day)
ds_t = open_any(row["tas_filepath"]).chunk({"lat": -1, "lon": -1, "time": 30})
ds_t = ds_t.sel(time = ds_t.time.dt.year.isin(years_insample))

pvar, ds_p = None, None
has_precip = ("precip_filepath" in row and isinstance(row["precip_filepath"], str)
              and row["precip_filepath"] and str(row["precip_filepath"]).lower() not in ("nan","none"))
if has_precip:
    ds_p = open_any(row["precip_filepath"]).chunk({"lat": -1, "lon": -1, "time": 30})
    ds_p = ds_p.sel(time = ds_p.time.dt.year.isin(years_insample))
    # MERRA2 prcorr→pr
    if product.upper() == "MERRA2" and ("prcorr" in ds_p.data_vars) and ("pr" not in ds_p.data_vars):
        ds_p = ds_p.rename({"prcorr": "pr"})
    pvar = detect_precip_var(ds_p)
    if pvar is not None:
        ds_p = normalize_precip_to_mm_per_day(ds_p, pvar, product)

ds_vars = [ds_t[["tas"]]]
if pvar is not None:
    ds_vars.append(ds_p[[pvar]])
ds = xr.merge(ds_vars).chunk({"time": 30})
ds = ds.assign(tas = ds.tas - 273.15)  # Kelvin → °C
ds["tas"].attrs["units"] = "degC"

# Weight maps
sample = ds.isel(time=0, drop=True)
with xa.set_options(silent=True):
    wm_adm2 = xa.pixel_overlaps(sample, gdf_sub)
    wm_adm1 = xa.pixel_overlaps(sample, gdf_adm1)

# Long-run T (ADM1): mean of yearly means, 1981–2010 (fallback to panel years)
clim_start, clim_end = 1981, 2010
years_lr = [y for y in years_insample if (clim_start <= y <= clim_end)] or list(years_insample)

ds_lr = ds.sel(time = ds.time.dt.year.isin(years_lr))
t_sum_yr  = ds_lr.tas.resample(time="1YS").sum(keep_attrs=True)
n_days_yr = ds_lr.tas.resample(time="1YS").count()
tmean_yr  = t_sum_yr / n_days_yr
tmean_lr  = tmean_yr.mean("time").chunk({"lat": -1, "lon": -1}).rename("Tmean")

with xa.set_options(impl="numba", silent=True):
    agg_lr = xa.aggregate(tmean_lr, wm_adm1)
df_lr = (
    agg_lr.to_dataframe().reset_index()[["iso","adm1_id","Tmean"]]
         .rename(columns={"Tmean": f"lr_tavg_{TAG}_adm1_avg"})
)

# Yearly aggregation (sum of daily powers)
rows = []
for yr in years_insample:
    print(f"[INFO] Aggregating year {yr}")
    ds_y = ds.sel(time = ds.time.dt.year == yr)

    t1 = ds_y.tas.sum("time").rename("T1_sum")
    t2 = (ds_y.tas**2).sum("time").rename("T2_sum")
    t3 = (ds_y.tas**3).sum("time").rename("T3_sum")
    t4 = (ds_y.tas**4).sum("time").rename("T4_sum")
    vars_to_merge = [t1, t2, t3, t4]

    if pvar is not None:
        pr1 = ds_y[pvar].sum("time").rename("PR1_sum")
        pr2 = (ds_y[pvar]**2).sum("time").rename("PR2_sum")
        pr3 = (ds_y[pvar]**3).sum("time").rename("PR3_sum")
        pr4 = (ds_y[pvar]**4).sum("time").rename("PR4_sum")
        vars_to_merge += [pr1, pr2, pr3, pr4]

    ds_poly = xr.merge(vars_to_merge).chunk({"lat": -1, "lon": -1})

    with xa.set_options(impl="numba", silent=True):
        agg_adm2_year = xa.aggregate(ds_poly, wm_adm2)

    df = (
        agg_adm2_year.to_dataframe().reset_index()
        .rename(columns={
            "T1_sum": f"tavg_poly_1_{TAG}",
            "T2_sum": f"tavg_poly_2_{TAG}",
            "T3_sum": f"tavg_poly_3_{TAG}",
            "T4_sum": f"tavg_poly_4_{TAG}",
        })
        .assign(year=int(yr), product=product)
    )

    if pvar is not None:
        df = df.rename(columns={
            "PR1_sum": f"prcp_poly_1_{TAG}",
            "PR2_sum": f"prcp_poly_2_{TAG}",
            "PR3_sum": f"prcp_poly_3_{TAG}",
            "PR4_sum": f"prcp_poly_4_{TAG}",
        })
    else:
        for k in (1,2,3,4):
            df[f"prcp_poly_{k}_{TAG}"] = 0.0

    df = df.merge(df_lr, on=["iso","adm1_id"], how="left")
    rows.append(df)

out = pd.concat(rows, ignore_index=True)

# Final EU collapse (ids already panel-aligned)
out["iso"] = collapse_eu_iso(out["iso"])

# PANEL-ANCHORED: guarantee every panel row exists
universe = keys[["iso","adm1_id","adm2_id","year"]].drop_duplicates()
out = universe.merge(out, on=["iso","adm1_id","adm2_id","year"], how="left")

# Tidy dtypes
out["adm1_id"] = out["adm1_id"].astype(str).str.strip()
out["adm2_id"] = out["adm2_id"].astype(str).str.strip()
out["year"]    = out["year"].astype(int)

out.to_csv(out_csv, index=False)
print(f"✅ Wrote {out_csv.name} with {len(out):,} rows at {OUTDIR}")
print("[SANITY] Years:", int(out['year'].min()), "…", int(out['year'].max()))
print("[SANITY] Non-null tavg_poly_1:", out[f'tavg_poly_1_{TAG}'].notna().sum())
print("[SANITY] IND rows:", len(out.loc[out['iso']=='IND']))
