#!/bin/bash
#
# Generate nudging files for CRITFC CMOP SELFE forecasts from liveocean data - uses anaconda python 2.7 
# library dependencies are:
# netCDF4, numpy, scipy, 
# dateutil, argparse, os, sys

unset DISPLAY

export ANACONDA=/usr/local/anaconda
export PYTHONPATH=$ANACONDA/lib/python2.7/site-packages
export LD_LIBRARY_PATH=$ANACONDA/lib:$LD_LIBRARY_PATH
export PATH=$ANACONDA/bin:$PATH
export BINDIR=$ANACONDA/bin

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

