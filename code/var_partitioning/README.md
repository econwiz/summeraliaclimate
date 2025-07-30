# var_partitioning
 **Decomposes variance** for mortality and GDP impacts across different GWL scenarios into **model** vs **internal** components, and to **generate summary figures**.

## Files

- **data.py**  
  Data‐loading and preprocessing utilities.  
  Reads your ensemble outputs (NetCDF/CSV) and prepares them for variance analysis.

- **stats.py**  
  Core variance‐partitioning routines.  
  - `calculate_impact_test(gwl)`: load or run the ensemble at a given global‑warming level.  
  - `variability_partitioning(ds)`: split each grid‑cell’s total variance into “model” vs “internal.”

- **plot.py**  
  Plotting helpers built on Matplotlib.  
  - Functions to draw bar charts, absolute‐variance panels, and other standard layouts.

- **gen_abs_var_figs.py**  
  Stand‑alone script that produces **absolute variance** figures (in original units) for each impact variable.

- **gen_bar_graphs.py**  
  Script to produce **fractional** bar charts showing the share of total variance explained by model vs internal variability.

- **generate_figures.py**  
  Master driver: runs the stats and plot modules in sequence to regenerate **all** figures in one go.

## Quick Start

1. Activate your environment:
   ```bash
   conda activate hle_iv


