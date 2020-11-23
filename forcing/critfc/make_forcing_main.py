"""
This is the main program for making the CRITFC output file,
using code from Charles Seaton.

For testing on my mac run in ipython as
run make_forcing_main.py -d 2019.07.04

"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

import subprocess
from datetime import datetime
start_time = datetime.now()

out_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + Ldir['date_string'] + '/'
out_fn = out_dir + 'critf.nc'

# critfc code inputs
rundate = Ldir['date_string'].replace('.','-')
depthfile = out_dir + 'ocean_his_0001.nc'
hgrid = Ldir['data'] + 'critfc/hgrid.ll'
vgrid = './vgrid.in'
basedir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
outdir = out_dir

    
cmd = ['python', 'gen_cmop_nudge.py', hgrid, vgrid, depthfile, basedir, outdir, rundate]
proc = subprocess.Popen(cmd)
proc.communicate()

"""
# run starts on today
rundate=`date --date "today" +%Y-%m-%d`
# file containing lon_rho, lat_rho, and depths variables
depthfile=/home/workspace/ccalmr46/liveocean/ocean_depths_20190712.nc
# CRITFC SELFE horizontal grid file
hgrid=/home/corie/forecasts/f33wrf/today/run/hgrid.ll
# CRITFC SELFE vertical grid file
vgrid=/home/corie/forecasts/f33wrf/today/run/vgrid.in
# base directory for LiveOcean files, assumes data files are in daily files below this base directory
basedir=/home/workspace/ccalmr46/liveocean/
# directory to write output CRITFC SELFE nudging files
outdir=/tmp

# python script path
scriptdir=./

cd $scriptdir

pwd
echo $BINDIR/python gen_cmop_nudge.py $hgrid $vgrid $depthfile $basedir $outdir $rundate
$BINDIR/python gen_cmop_nudge.py $hgrid $vgrid $depthfile $basedir $outdir $rundate

"""


#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

#ffun.finale(result_dict, Ldir, Lfun)