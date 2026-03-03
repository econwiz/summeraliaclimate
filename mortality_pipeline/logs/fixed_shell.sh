#!/user/ab5405/.conda/envs/hle_iv/bin/python
import subprocess
from dict_products import products

for product in products:
    subprocess.run([
      "grid_run",
      "--grid_submit=batch",
      "--grid_mem=16G",
      "--grid_email=ab5405@columbia.edu",
      "../code/fixed_agg_to_run.py",
      "--product", product
    ])
