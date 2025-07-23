import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pathlib import Path
import cartopy.crs as ccrs
import xarray as xr

#simple map plotting function, plots CONUS with borders, using the reds colorscale
def plot_map(data, title, cbar_label):
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    data.plot(ax=ax, cmap='Reds', cbar_kwargs={'label': cbar_label})
    ax.coastlines()
    ax.set_title(title)
    plt.show()

# Create a two‐row panel of maps showing how internal and model variability contribute to local climate variability across different global warming levels (GWLs). Each subplot displays the fractional contribution for one source (internal or model) at a given GWL, with consistent color scaling and labeled panels. If an output directory is provided, the figure is saved; otherwise it is displayed inline.
def plot_variability_panels(stats_ds: xr.Dataset,
                            title: str = "Variability partitioning",
                            cmap: str = 'Greens',
                            cbar_label: str = 'Fractional contribution to local variability',
                            outdir: Path = None,
                            fname_prefix: str = 'panel') -> None:
    
    gwls    = stats_ds.gwl.values
    sources = ['internal', 'model']
    ncols   = len(gwls)
    
    fig, axes = plt.subplots(nrows=2, ncols=ncols,figsize=(2*ncols, 8), subplot_kw={'projection': ccrs.PlateCarree()},)
    plt.subplots_adjust(hspace=0.3, wspace=0.1)
    
    panel_labels = [f"{chr(97 + k)}." for k in range(2 * ncols)]
    
    for i, src in enumerate(sources):
        for j, gwl in enumerate(gwls):
            ax = axes[i, j]
            da = stats_ds[f"f_{src}"].sel(gwl=gwl)
            
            im = da.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(),cmap=cmap, vmin=0.0, vmax=1.0,add_colorbar=False)
            ax.coastlines(color='gray', linewidth=0.5)
            ax.set_xticks([]); ax.set_yticks([])
            
            if i == 0:
                ax.set_title(f"{gwl:.1f}°C", fontsize=12)
            if j == 0:
                ax.text(-0.12, 0.5, src.capitalize(),
                    transform=ax.transAxes,
                    va='center', ha='right',
                    rotation=90, fontsize=12
                )
                
            ax.text(
                0.02, 0.98, panel_labels[i * ncols + j],
                transform=ax.transAxes,
                va='top', ha='left',
                fontweight='bold')
    
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])

    sm = plt.cm.ScalarMappable(
        norm=plt.Normalize(vmin=0.0, vmax=1.0),
        cmap=cmap
    )
    sm.set_array([])
    fig.colorbar(sm, cax=cbar_ax, label=cbar_label)
    
    fig.suptitle(title, fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    
    if outdir is not None:
        filepath = outdir / f"{fname_prefix}_{ncols}gwl_panels.png"
        fig.savefig(filepath, dpi=100, bbox_inches='tight')
        plt.close(fig)   # <- prevent inline display
    else:
        plt.show()
        
        