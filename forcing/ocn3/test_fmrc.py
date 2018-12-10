""""
Code to test getting hycom files using the new FMRC_best file.
"""

import netCDF4 as nc

fn = ('http://tds.hycom.org/thredds/dodsC/' +
        'GLBu0.08/expt_93.0/data/forecasts/FMRC_best.ncd')
        
ds = nc.Dataset(fn)

tvec = ds['time'][:]

for t in tvec[:3]:
    print(t)
    
ds.close()