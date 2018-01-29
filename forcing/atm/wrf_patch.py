"""
Code to replace missing WRF files with those from another day.

This will just copy all the files from the substitute day to
the missing day, changing their names to be consistent with
the missing day.  This does not change their internal times,
but the cludgey way I make ROMS atm forcing files creates its
own time variable so it does not matter.

"""

import os
import shutil

dir0="/pmr2/darr/wrf_crons/wrfout/"

substitute_day = "2017092200"

missing_day =    "2017092300"

indir = dir0 + substitute_day + '/'
outdir = dir0 + missing_day + '/'

# I assume that outdir exists, is empty,
# and that I have write permission there

in_list = os.listdir(indir)

for in_fn in in_list:
    if in_fn[:6] == 'wrfout':
        
        out_fn = in_fn.replace(substitute_day, missing_day)
        
        if False: # Testing
            print('Copying ' + indir+in_fn)
            print('-- to - ' + outdir+out_fn)
            
        else:
            shutil.copyfile(indir+in_fn, outdir+out_fn)


