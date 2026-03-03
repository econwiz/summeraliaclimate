#!/bin/bash
module load anaconda/2023.07
source activate hle_iv

cd /home/ab5405/summeraliaclimate/code
./agg_to_reg.py
