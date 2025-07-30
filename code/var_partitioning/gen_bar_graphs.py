#!/usr/bin/env python3
import sys
from pathlib import Path

code_dir = Path(__file__).resolve().parent
if str(code_dir) not in sys.path:
    sys.path.insert(0, str(code_dir))

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#calculate_impact_test returns ensemble arrays for 1 GWL var_partitioning decomposes variance into model and internal
from var_partitioning.stats import calculate_impact_test, variability_partitioning

script_path = Path(__file__).resolve()
code_dir    = script_path.parent 


def main():
    #get the 5 gwls
    gwl_list = [1.0, 1.5, 2.0, 2.5, 3.0]
    
    #Split ensemble variance into model vs internal variability using the variability_partitioning function. First 
    #calculate the mortality and gdp impact arrays for each GWL. 
    mort_list, gdp_list = [], []
    for g in gwl_list:
        mort_all, gdp_all = calculate_impact_test(g)
        mort_ds = variability_partitioning(mort_all).assign_coords(gwl=g)
        gdp_ds  = variability_partitioning(gdp_all).assign_coords(gwl=g)
        mort_list.append(mort_ds)
        gdp_list.append(gdp_ds)

    #Concatenate the list of per‑GWL mortality datasets along a new 'gwl' dimension
    mort_stats = xr.concat(mort_list, dim="gwl").load()
    gdp_stats  = xr.concat(gdp_list,  dim="gwl").load()

    #Create a boolean mask to identify grid cells where the total variance is positive
    mask_mort = mort_stats["total_variance"] > 0
    mask_gdp = gdp_stats["total_variance"] > 0 

    #Select the 'f_internal' fraction only on cells with positive total variance from the var_partitioning function
    #compute the spatial mean over lat & lon, extract result as a numpy array
    f_int_mort = mort_stats["f_internal"].where(mask_mort).mean(dim=("lat","lon"), skipna=True).values
    f_mod_mort = mort_stats["f_model"].where(mask_mort).mean(dim=("lat","lon"), skipna=True).values
    f_int_gdp  = gdp_stats["f_internal"].where(mask_gdp).mean(dim=("lat","lon"), skipna=True).values
    f_mod_gdp  = gdp_stats["f_model"].where(mask_gdp).mean(dim=("lat","lon"), skipna=True).values

    figures_dir = code_dir.parent / "figures"

    #Plot stacked‐bar chart of mortality variance fractions (model vs internal) across GWLs for mortality and gdp. 
    #Save to a specific figures directory
    fig, ax = plt.subplots(figsize=(6,4))
    ax.bar(gwl_list, f_mod_mort, width=0.4, label="model")
    ax.bar(gwl_list, f_int_mort, width=0.4, bottom=f_mod_mort, label="internal")
    ax.set_xlabel("Global-warming level (°C)")
    ax.set_ylabel("Fraction of variance")
    ax.set_title("Mortality uncertainty fractions by GWL")
    ax.set_xticks(gwl_list)
    ax.legend()
    plt.tight_layout()
    if not figures_dir.exists():
        figures_dir.mkdir(parents=True)
    fig.savefig(figures_dir / "mortality_fraction_bar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6,4))
    ax.bar(gwl_list, f_mod_gdp, width=0.4, label="model")
    ax.bar(gwl_list, f_int_gdp, width=0.4, bottom=f_mod_gdp, label="internal")
    ax.set_xlabel("Global-warming level (°C)")
    ax.set_ylabel("Fraction of variance")
    ax.set_title("GDP uncertainty fractions by GWL")
    ax.set_xticks(gwl_list)
    ax.legend()
    plt.tight_layout()
    if not figures_dir.exists():
        figures_dir.mkdir(parents=True)
    fig.savefig(figures_dir / "gdp_fraction_bar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    
    print("Bar charts saved to", figures_dir)

if __name__ == "__main__":
    main()
