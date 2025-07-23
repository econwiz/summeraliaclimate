import pandas as pd
import numpy as np
import xarray as xr

from .data import rebin_temp_data_mortality, rebin_temp_data_gdp
from funcs_support import get_filepaths

#Compute the partitioning of total variability into internal variability and model uncertainty for a multi-model ensemble of impact projections. Internal variability is calculated as the run‐to‐run variance within each model and averaged across models, while model uncertainty is the variance of model‐mean impacts across models. The function returns a Dataset containing the internal variance, model variance, total variance, and their fractional contributions for each grid cell or region.
def variability_partitioning(impact_all):
    impact_all = impact_all.where(impact_all != 0)
    
    #this finds the internal variability within each model over different runs
    within_model_var = impact_all.groupby('model').var(dim='run', skipna=True, ddof=0)
    #finds mean internal variability over the dimension model (i.e. mean internal var across models)
    internal_var = within_model_var.mean(dim='model', skipna=True)

    #this finds the mean impacts by model, for each run 
    model_means = impact_all.groupby('model').mean(dim='run', skipna=True)
    #this finds the variance in model means across models (i.e. model uncertainty)
    model_unc = model_means.var(dim='model', skipna=True, ddof=0)

    #calculates total variability as the sum of internal and model uncertainty
    total_var = internal_var + model_unc

    #gets fractional internal and model uncertainty, making sure to avoid dividing by zeros
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

#calculates the impact (either in terms of mortality or GDP percent change) for each additional day 
#spent in a given temperature bin at the grid cell level, with gwl as the input. Thus, calculates a
#map of impacts over CONUS at the grid-cell level.
def calculate_impact_test(gwl):
    # Prepare lists to collect per-run outputs
    mort_impact_list = []
    gdp_impact_list  = []
    model_labels     = []
    run_labels       = []
    
    #iterates through each entry and row in the dataframe containing temperature data, and finds the 
    #change in the number of days at GWL 0.61 and the input gwl. The impact functions are rebinned
    #according to the temperature data, and the difference in days in each bin is multiplied by the 
    #imapct coefficeint. These impacts are then added to the empty lists created above, and concatenated
    #on a run dimension.
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