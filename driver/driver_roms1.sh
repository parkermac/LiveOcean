#!/bin/bash

# This runs ROMS for one or more days, allowing for either
# a forecast or backfill.
#
# CHANGES:
# 8/2/2016 This version is designed to handle ROMS jobs more gracefully,
# using advice from here:
# http://unix.stackexchange.com/questions/76717/bash-launch-background-process-and-check-when-it-ends

# NOTE: must be run from gaggle, and depends on other drivers having been run first

# set paths and connect to a library of functions
if [ $HOME = "/Users/PM5" ] ; then
  LO_parent="/Users/PM5/Documents/LiveOcean"
  R_parent="/Users/PM5/Documents/LiveOcean_roms"
elif [ $HOME = "/home/parker" ] ; then
  LO_parent="/fjdata1/parker/LiveOcean"
  R_parent="/pmr1/parker/LiveOcean_roms"
fi
. $LO_parent"/driver/common.lib"

# USE COMMAND LINE OPTIONS
#
# -g name of the grid [cascadia1, ...]
# -t name of the forcing tag [base, ...]
# -x name of the ROMS executable to use
# -s new or continuation
# -r forecast or backfill
#  if backfill, then you must provide two more arguments:
# -0 start date: yyyymmdd
# -1 end date: yyyymmdd
#
# example call to do backfill, with a cold start:
# ./driver_roms1.sh -g cascadia1 -t base -x lobio1 -s new -r backfill -0 20140201 -1 20140203
#
# example call to do forecast:
# ./driver_roms1.sh -g cascadia1 -t base -x lobio1 -s continuation -r forecast
#
# you can also use long names like --ex_name instead of -x

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
    -s | --start_type )  shift
      start_type=$1
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

# initialize control flags
keep_going=1 # 1 => keep going, 0 => stop the driver
blow_ups=0

# start the main loop over days
while [ $D -le $D1 ] && [ $keep_going -eq 1 ]
do
  echo "********** driver_roms1.sh *******************"
  echo "  blow ups = " $blow_ups
  
  # manipulate the string D to insert dots, using the syntax:
  # substring = ${string:start_index:count} and the index starts from 0
  DD=${D:0:4}.${D:4:2}.${D:6:2}

  f_string="f"$DD
  Rf=$R_parent"/output/"$gtagex"/"$f_string
  log_file=$Rf"/log"

  # choose which cores to run on (must be consistent in dot_in)
  if [ $run_type = "forecast" ] ; then
    hf="../shared/hf144"
    np_num=144
  elif [ $run_type = "backfill" ] ; then
    hf="../shared/hf144"
    np_num=144
  fi

  # Run make_dot_in.py, which creates an empty f_string directory,
  # and then cd to where the ROMS executable lives.
  cd $LO_parent"/forcing/dot_in/"$gtagex
  source $HOME"/.bashrc"
  if [ $D = $D0 ] && [ $start_type = "new" ] ; then
    python ./make_dot_in.py -g $gridname -t $tag -s $start_type -r $run_type -d $DD -x $ex_name -np $np_num -bu $blow_ups
    sleep 30
    cd $R_parent"/makefiles/"$ex_name"_tideramp"
  else
    python ./make_dot_in.py -g $gridname -t $tag -s continuation -r $run_type -d $DD -x $ex_name -np $np_num -bu $blow_ups
    sleep 30
    cd $R_parent"/makefiles/"$ex_name
  fi

  # the actual ROMS run command
  if [ $HOME = "/Users/PM5" ] ; then # testing
    echo "/cm/shared/local/openmpi-ifort/bin/mpirun -np $np_num -machinefile $hf oceanM $Rf/liveocean.in > $log_file &"
  elif [ $HOME == "/home/parker" ] ; then # the real thing
    /cm/shared/local/openmpi-ifort/bin/mpirun -np $np_num -machinefile $hf oceanM $Rf/liveocean.in > $log_file &
    # Check that ROMS has finished successfully.
    PID1=$!
    wait $PID1
    echo "run completed for" $f_string
  fi

  # check the log_file to see if we should continue
  if grep -q "Blowing-up" $log_file ; then
    echo "- Run blew up!"
    blow_ups=$(( $blow_ups + 1 )) #increment the blow ups
    if [ $blow_ups -le 1 ] ; then
      keep_going=1
    else
      keep_going=0
    fi
  elif grep -q "ERROR" $log_file ; then
    echo "- Run had an error."
    keep_going=0
  elif grep -q "ROMS/TOMS: DONE" $log_file ; then
    echo "- Run completed successfully."
    keep_going=1
    blow_ups=0
  else
    echo "- Something else happened."
    keep_going=0
  fi

  echo $(date)

  if [ $blow_ups -eq 0 ] ; then
    # This function increments the day.
    # NOTE: it changes y, m, d, and D, even in the scope of this shell script!
    next_date $y $m $d
    blow_ups=0 # reset the blow up counter
  fi

done # end of while loop

