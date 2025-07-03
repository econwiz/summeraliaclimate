import pandas as pd
import numpy as np
import xarray as xr

from .data import rebin_temp_data_mortality, rebin_temp_data_gdp
from funcs_support import get_filepaths

def variability_partitioning(impact_all):
    impact_all = impact_all.where(impact_all != 0)
    
    within_model_var = impact_all.groupby('model').var(dim='run', skipna=True, ddof=0)
    internal_var = within_model_var.mean(dim='model', skipna=True)
    
    model_means = impact_all.groupby('model').mean(dim='run', skipna=True)
    model_unc = model_means.var(dim='model', skipna=True, ddof=0)
    
    total_var = internal_var + model_unc
    f_internal = xr.where(total_var > 0, internal_var / total_var, 0.0)
    f_model    = xr.where(total_var > 0, model_unc    / total_var, 0.0)
    
    return xr.Dataset({
        "internal_variance": internal_var,
        "model_unc":         model_unc,
        "total_variance":    total_var,
        "f_internal":        f_internal,
        "f_model":           f_model
    })

dir_list_df = pd.read_csv(
    '/user/ab5405/summeraliaclimate/code/dir_list.csv',
    index_col=0
)
dir_list = dir_list_df['dir_path'].to_dict()
df_paths = get_filepaths(source_dir='proc', dir_list=dir_list)

df = df_paths.query("varname == 'tasdmgfparams' and gwl == 'ALLGWLs'").copy()

def calculate_impact_test(gwl):
    # Prepare lists to collect per-run outputs
    mort_impact_list = []
    gdp_impact_list  = []
    model_labels     = []
    run_labels       = []

    for _, row in df.iterrows():
        model_labels.append(row['model'])
        run_labels.append(row['run'])
    
        ds = xr.open_zarr(row['path'], consolidated=False)
        
        gwl_delta = ds['bins_behrer'].interp(gwl=gwl).transpose('bin_f','lat','lon')
        gwl_06    = ds['bins_behrer'].interp(gwl=0.61).transpose('bin_f','lat','lon')
        diff = gwl_delta - gwl_06
    
        mort_coef_fine = rebin_temp_data_mortality(diff)
        mort_map = (diff * mort_coef_fine).sum(dim='bin_f')
        mort_impact_list.append(mort_map)

        gdp_coef_fine = rebin_temp_data_gdp(diff)
        gdp_map_ln = (diff * gdp_coef_fine).sum(dim='bin_f')
        gdp_all_pct = (np.exp(gdp_map_ln) - 1) * 100
        gdp_impact_list.append(gdp_all_pct)
        
    mort_all = xr.concat(mort_impact_list, dim='run')
    mort_all = mort_all.assign_coords({
        'run':   run_labels,
        'model': ('run', model_labels)
    })
    
    gdp_all = xr.concat(gdp_impact_list, dim='run')
    gdp_all = gdp_all.assign_coords({
        'run':   run_labels,
        'model': ('run', model_labels)
    })
    
    return mort_all, gdp_all