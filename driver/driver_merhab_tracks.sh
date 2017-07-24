#!/bin/bash

# This runs the code to create a figre of particle tracks for MERHAB.

# NOTE: must be run from fjord.

# set a path and connect to a library of functions
if [ $HOME = "/Users/PM5" ] ; then
  LO_parent="/Users/PM5/Documents/LiveOcean"
elif [ $HOME = "/home/parker" ] ; then
  LO_parent="/data1/parker/LiveOcean"
fi
. $LO_parent"/driver/common.lib"

# USE COMMAND LINE OPTIONS
#
# -g name of the grid [cascadia1, ...]
# -t name of the forcing tag [base, ...]
# -x name of the ROMS executable to use [lobio1, ...]
# -r run type [forecast, backfill]
#  if backfill, then you must provide
# -0 start date: yyyymmdd
#
# example call
#
# ./driver_merhab_tracks.sh -g cascadia1 -t base -x lobio1 -r forecast
# ./driver_merhab_tracks.sh -g cascadia1 -t base -x lobio1 -r backfill -0 20170518

ex_name="placeholder"
while [ "$1" != "" ] ; do
  case $1 in
    -g | --gridname )  shift
      gridname=$1
      ;;
    -t | --tag )  shift
      tag=$1
      ;;
    -x | --ex_name )  shift
      ex_name=$1
      ;;
    -r | --run_type )  shift
      run_type=$1
      ;;
    -0 | --ymd0 )  shift
      ymd0=$1
      ;;
  esac
  shift
done

if [ $run_type = "forecast" ] ; then
  # do forecast
  ymd0=$(date "+%Y%m%d")
fi

# parse the date strings
y0=${ymd0:0:4}; m0=${ymd0:4:2}; d0=${ymd0:6:2}

y=$y0; m=$m0; d=$d0
# note the leading "10#" so that it doesn't interpret 08 etc. as octal
D0=$[10#$y*10000 + 10#$m*100 + 10#$d]
D=$D0

# manipulate the string D to insert dots, using the syntax:
#    substring = ${string:start_index:count}
# and the index starts from 0
DD=${D:0:4}.${D:4:2}.${D:6:2}

# Make the forcing.
cd $LO_parent"/plotting/"
source $HOME"/.bashrc"
if [ -e $HOME"/.bash_profile" ] ; then
  source $HOME"/.bash_profile"
fi
if [ -e $HOME"/.profile" ] ; then
  source $HOME"/.profile"
fi

# could add -mov True to this
python ./pan_plot.py -g $gridname -t $tag -x $ex_name -lt merhab -d $DD -pt P_tracks_MERHAB -mov True &

