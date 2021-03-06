#!/bin/bash

# Code to test the date increment and decrement functions in common.lib

# Required command line arguments
#
# -0 start date: yyyymmdd
# -1 end date: yyyymmdd
#
# Example call
# ./test_date.sh -0 20200226 -1 20200304

# run the code to put the environment into a csv
if [ -e ../../LiveOcean_user/alpha/get_lo_info.sh ] ; then
  ../../LiveOcean_user/alpha/get_lo_info.sh
else
  ../alpha/get_lo_info.sh
fi

# and read the csv into active variables
while IFS=, read col1 col2
do
  eval $col1="$col2"
done < ../alpha/lo_info.csv

. $LO"driver/common.lib"

while [ "$1" != "" ]; do
  case $1 in
    -0 | --ymd0 )  shift
      ymd0=$1
      ;;
    -1 | --ymd1 )  shift
      ymd1=$1
      ;;
  esac
  shift
done

# parse the date strings
y0=${ymd0:0:4}; m0=${ymd0:4:2}; d0=${ymd0:6:2}
y1=${ymd1:0:4}; m1=${ymd1:4:2}; d1=${ymd1:6:2}

y=$y0; m=$m0; d=$d0
# note the leading "10#" so that it doesn't interpret 08 etc. as octal
D0=$[10#$y*10000 + 10#$m*100 + 10#$d]
D=$D0
D1=$[10#$y1*10000 + 10#$m1*100 + 10#$d1]

# start the main loop over days
while [ $D -le $D1 ]
do
  # This function changes the value in $DPP to be two days before that in $D
  two_days_before $y $m $d
  echo "Day = "$D" (Two Days Before = "$DPP")"

  # This function increments the day.
  # NOTE: it changes y, m, d, and D, even in the scope of this shell script!
  next_date $y $m $d
done # end of while loop