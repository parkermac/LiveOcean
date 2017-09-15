#!/bin/bash

# Code to run all the post-processing jobs for a forecast
# must be run drom LiveOcean/driver/

# set a path and connect to a library of functions
if [ $HOME = "/Users/pm7" ] ; then
  LO_parent=$HOME"/Documents/LiveOcean"
elif [ $HOME = "/home/parker" ] ; then
  LO_parent="/data1/parker/LiveOcean"
fi

for frc in 'tracks_m' 'carbon' 'low_pass' 'ubc' 'surface' 'azu1'; do

  ./driver_forcing1.sh -g cascadia1 -t base -x lobio1 -f $frc -r forecast > $LO_parent/driver/dlog_$frc &
  # Check that the job has finished successfully.
  PID1=$!
  wait $PID1
  echo "job completed for" $frc "at" $(date)
  echo $(date)
  sleep 10

done