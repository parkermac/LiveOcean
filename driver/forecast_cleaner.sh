#!/bin/bash

# This deletes the day 2 and day 3 history files from past forecasts.
# Meant to be run by cron late every day, well after the files are needed.
# Hard-wired only to work on boiler, with selected runs.

# Also deletes some other post-processing products.

rm /data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f*/ocean_his_002[6-9].nc
rm /data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f*/ocean_his_00[3-7][0-9].nc

rm /data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f*/cmop_*.nc
rm /data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f*/ocean_layers.nc
rm /data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f*/ocean_surface.nc
rm /data1/parker/LiveOcean_roms/output/cas6_v3_lo8b/f*/low_passed_UBC.nc
