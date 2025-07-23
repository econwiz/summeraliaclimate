#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path

from var_partitioning.stats import calculate_impact_test
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def plot_mean_impacts(mort_means, gdp_means, gwl_list, outdir):
    n = len(gwl_list)
    # make each subplot ~3" wide and total height ~6"
    fig, axes = plt.subplots(
        2, n, figsize=(3*n, 6),
        subplot_kw={"projection": ccrs.PlateCarree()},
        gridspec_kw={"wspace": 0.1, "hspace": 0.2, "right": 0.88}
    )

    # determine vmin/vmax for each row
    mort_vmin = mort_means.min().item()
    mort_vmax = mort_means.max().item()
    gdp_vmin  = gdp_means .min().item()
    gdp_vmax  = gdp_means .max().item()

    # plot mortality (row 0)
    for i, gwl in enumerate(gwl_list):
        ax = axes[0, i]
        im1 = mort_means.sel(gwl=gwl).plot(
            ax=ax, cmap="Reds", vmin=mort_vmin, vmax=mort_vmax,
            add_colorbar=False
        )
        ax.coastlines()
        ax.set_extent([-125, -66.5, 24, 50], ccrs.PlateCarree())
        ax.set_title(f"{gwl}°C")
        if i == 0:
            ax.set_ylabel("Mortality")

    # plot GDP (row 1)
    for i, gwl in enumerate(gwl_list):
        ax = axes[1, i]
        im2 = gdp_means.sel(gwl=gwl).plot(
            ax=ax, cmap="Purples", vmin=gdp_vmin, vmax=gdp_vmax,
            add_colorbar=False
        )
        ax.coastlines()
        ax.set_extent([-125, -66.5, 24, 50], ccrs.PlateCarree())
        if i == 0:
            ax.set_ylabel("GDP per capita")

    # add colorbars in reserved space on the right
    cbar_ax1 = fig.add_axes([0.90, 0.56, 0.02, 0.33])  # [left, bottom, width, height]
    fig.colorbar(im1, cax=cbar_ax1, label="Mean mortality impact")

    cbar_ax2 = fig.add_axes([0.90, 0.11, 0.02, 0.33])
    fig.colorbar(im2, cax=cbar_ax2, label="Mean GDP impact")

    plt.suptitle("Mean mortality and GDPpc impacts across GWLs", y=0.95)
    plt.tight_layout(rect=[0, 0, 0.88, 0.93])

    outdir.mkdir(parents=True, exist_ok=True)
    fpath = outdir / "mean_impacts_gwl_panel_fixed2.png"
    plt.savefig(fpath, dpi=300, bbox_inches="tight")
    print(f"Saved fixed figure: {fpath}")

def main():
    # locate code & figures folder
    script_path = Path(__file__).resolve()
    code_dir    = script_path.parent
    figures_dir = code_dir.parent / "figures"

    gwl_list = [1.0, 1.5, 2.0, 2.5, 3.0]

    # collect mean-impact arrays
    mort_means = []
    gdp_means  = []
    for gwl in gwl_list:
        print(f"→ Processing GWL = {gwl}°C …")
        mort_all, gdp_all = calculate_impact_test(gwl)

        mort_means.append(
            mort_all.mean(dim="run").assign_coords(gwl=gwl)
        )
        gdp_means. append(
            gdp_all.mean(dim="run").assign_coords(gwl=gwl)
        )

    # concatenate along the GWL dimension and load into memory
    mort_means = xr.concat(mort_means, dim="gwl").load()
    gdp_means  = xr.concat(gdp_means,  dim="gwl").load()

    # plot & save
    plot_mean_impacts(mort_means, gdp_means, gwl_list, figures_dir)

if __name__ == "__main__":
    main()
