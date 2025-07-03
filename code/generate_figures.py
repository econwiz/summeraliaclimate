#!/usr/bin/env python3
import sys, traceback
from pathlib import Path

def main():
    print("Importing your functions…")
    from var_partitioning.stats import calculate_impact_test, variability_partitioning
    from var_partitioning.plot  import plot_variability_panels

    import xarray as xr

    gwl_list = [1.0, 1.5, 2.0, 2.5, 3.0]
    print("GWL list:", gwl_list)
    
    mort_list, gdp_list = [], []
    for g in gwl_list:
        print(f"  → Processing GWL={g}°C …")
        mort_all, gdp_all = calculate_impact_test(g)
        mort_list.append(variability_partitioning(mort_all).assign_coords(gwl=g))
        gdp_list .append(variability_partitioning(gdp_all).assign_coords(gwl=g))

    print("Concatenating and loading into memory…")
    mort_stats = xr.concat(mort_list, dim="gwl").load()
    gdp_stats  = xr.concat(gdp_list,  dim="gwl").load()

    figures_dir = code_dir.parent / "figures"
    print(f"Figures will be saved to: {figures_dir}")
    if not figures_dir.exists():
        print("  (But that folder didn’t exist—creating it now.)")
        figures_dir.mkdir(parents=True)

    print("Saving mortality panel…")
    plot_variability_panels(
        mort_stats,
        title="County-level mortality variability partitioning",
        cmap="Greens",
        cbar_label="Fractional contribution to local variability",
        outdir=figures_dir,
        fname_prefix="mortality"
    )

    print("Saving GDP panel…")
    plot_variability_panels(
        gdp_stats,
        title="County-level per-capita GDP variability partitioning",
        cmap="Blues",
        cbar_label="Fractional contribution to local variability",
        outdir=figures_dir,
        fname_prefix="gdp"
    )

    print("All panels saved!")

if __name__ == "__main__":
    main()
