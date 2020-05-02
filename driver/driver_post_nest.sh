#!/bin/bash

# This runs post-processing jobs for a forecast, with a while loop
# that makes sure all of the history files are in place before starting.

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


# USE COMMAND LINE OPTIONS
#
# -g name of the grid [cascadia1, ...]
# -t name of the forcing tag [base, ...]
# -x name of the ROMS executable to use
# -r forecast or backfill
#  if backfill, then you must provide two more arguments:
# -0 start date: yyyymmdd
# -1 end date: yyyymmdd
#
# example call to do backfill
# NOTE: I don't think this will be called for backfill.  It is really
# meant to be used in forecast mode, called by cron.
#
# example call to do forecast:
# ./driver_post.sh -g cas4 -t v2 -x lo6biom -r forecast

while [ "$1" != "" ]; do
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
    -1 | --ymd1 )  shift
      ymd1=$1
      ;;
  esac
  shift
done

if [ $run_type = "forecast" ]
then
  # do forecast
  ymd0=$(date "+%Y%m%d")
  ymd1=$ymd0
fi

# parse the date strings
y0=${ymd0:0:4}; m0=${ymd0:4:2}; d0=${ymd0:6:2}
y1=${ymd1:0:4}; m1=${ymd1:4:2}; d1=${ymd1:6:2}

y=$y0; m=$m0; d=$d0
# note the leading "10#" so that it doesn't interpret 08 etc. as octal
D0=$[10#$y*10000 + 10#$m*100 + 10#$d]
D=$D0
D1=$[10#$y1*10000 + 10#$m1*100 + 10#$d1]

gtag=$gridname"_"$tag
gtagex=$gtag"_"$ex_name

# start the main loop over days
while [ $D -le $D1 ]
do
  echo "********** driver_post_nest.sh *******************"
  
  # manipulate the string D to insert dots, using the syntax:
  # substring = ${string:start_index:count} and the index starts from 0
  DD=${D:0:4}.${D:4:2}.${D:6:2}

  f_string="f"$DD
  Rf=$roms"output/"$gtagex"/"$f_string"/"
  
  if [ -e $HOME"/.bashrc" ] ; then
    source $HOME"/.bashrc"
  elif [ -e $HOME"/.bash_profile" ] ; then
    source $HOME"/.bash_profile"
  elif [ -e $HOME"/.profile" ] ; then
    source $HOME"/.profile"
  fi

  
  if [ $run_type == 'forecast' ] ; then
    nhis=73
  elif [ $run_type == 'backfill' ] ; then
    nhis=25
  fi

  if [ $lo_env == "pm_mac" ] ; then # testing
    maxcount=3
    sleeptime=1
  elif [ $lo_env == "pm_boiler" ] ; then
    maxcount=480
    sleeptime=60
  fi

  echo "- Looking in "$Rf
  all_files_here=0
  count=0
  while [ $all_files_here -eq 0 ] && [ $count -le $maxcount ]
  do
    # test to see if all files are here
    all_files_here=1
    for i in $(seq 1 $nhis); do
      his_num='000'$i
      his_num=${his_num: -4}
      fn=$Rf'ocean_his_'$his_num'.nc'
      if [ ! -f $fn ] ; then
        all_files_here=0
      fi
    done
    echo "-- count = "$count
    echo "-- all_files_here = "$all_files_here
    count=$[10#$count + 1]
    if [ $all_files_here -eq 0 ]; then
      sleep $sleeptime
    else
      sleep 11
    fi
  done
  
  # run post processing
  if [ $lo_env == "pm_mac" ] ; then # testing
    echo "TESTING"
  elif [ $lo_env == "pm_boiler" ] && [ $all_files_here -eq 1 ]; then
    for frc in 'active_forecast_nest' ; do
      # echo "Would be working on "$frc
      ./driver_forcing2.sh -g $gridname -t $tag -x $ex_name -f $frc -r $run_type > $LO"driver/dlog_"$frc &
      # Check that the job has finished successfully.
      PID1=$!
      wait $PID1
      echo "job completed for" $frc "at" $(date)
      echo $(date)
      sleep 15
    done

  fi
  
  echo $(date)

  # This function increments the day.
  # NOTE: it changes y, m, d, and D, even in the scope of this shell script!
  next_date $y $m $d

done # end of while loop