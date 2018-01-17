#!/bin/bash

# Code to run all the post-processing jobs for a forecast
# must be run drom LiveOcean/driver/

# run the code to put the environment into a csv
if [ -e ../alpha/user_get_lo_info.sh ] ; then
  ../alpha/user_get_lo_info.sh
else
  ../alpha/get_lo_info.sh
fi

# and read the csv into active variables
while IFS=, read col1 col2
do
  eval $col1="$col2"
done < ../alpha/lo_info.csv

for frc in 'tracks_m' 'carbon' 'low_pass' 'ubc' 'surface' 'azu1'; do

  ./driver_forcing2.sh -g cascadia1 -t base -x lobio5 -f $frc -r forecast > $LO"driver/dlog_"$frc &
  # Check that the job has finished successfully.
  PID1=$!
  wait $PID1
  echo "job completed for" $frc "at" $(date)
  echo $(date)
  sleep 10

done