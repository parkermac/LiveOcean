#!/bin/bash

# This runs ROMS for one or more days, allowing for either
# a forecast or backfill.

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

# COMMAND LINE INPUTS
#
# 1. new or continuation
# 2. forecast or backfill
#  if backfill, then you must provide two more arguments
#  3. start date: yyyymmdd
#  4. end date: yyyymmdd
#
# example call to do backfill:
# ./driver_roms.sh new backfill 20140201 20140203 > log_roms &
# example call from cron:
# ./driver_roms.sh continuation forecast

start_type=$1
run_type=$2

# check command line inputs
if [ $# -eq 2 ] ; then
  if [ $start_type = "new" ] || [ $start_type = "continuation" ] && [ $run_type = "forecast" ]
  then
    # do forecast
    ymd0=$(date "+%Y%m%d")
    ymd1=$ymd0
  else
    echo "Bad input arguments (two)"
    exit
  fi
elif [ $# -eq 4 ] ; then
  if [ $start_type = "new" ] || [ $start_type = "continuation" ] && [ $run_type = "backfill" ]
  then
    # do backfill
    ymd0=$3 # first date
    ymd1=$4 # last date
  else
    echo "Bad input arguments (four)"
    exit
  fi
else
  echo "Need 2 or 4 input arguments"
  exit
fi

# parse the date strings
y0=${ymd0:0:4}; m0=${ymd0:4:2}; d0=${ymd0:6:2}
y1=${ymd1:0:4}; m1=${ymd1:4:2}; d1=${ymd1:6:2}

y=$y0; m=$m0; d=$d0
# note the leading "10#" so that it doesn't interpret 08 etc. as octal
D0=$[10#$y*10000 + 10#$m*100 + 10#$d]
D=$D0
D1=$[10#$y1*10000 + 10#$m1*100 + 10#$d1]

# read gtag from the RUN_INFO.csv file
while read LINE
do
  KEY=${LINE%,*}
  VAL=${LINE##*,}
  if [ $KEY = "gridname" ] ; then
    gridname=$VAL
  elif [ $KEY = "tag" ] ; then
    tag=$VAL
  fi
done <$LO_parent"/alpha/RUN_INFO.csv"
gtag=$gridname"_"$tag

while [ $D -le $D1 ]
do
  # manipulate the string D to insert dots, using the syntax:
  #    substring = ${string:start_index:count}
  # and the index starts from 0
  DD=${D:0:4}.${D:4:2}.${D:6:2}

  f_string="f"$DD
  Rf=$R_parent"/output/"$gtag"/"$f_string
  #echo $Rf
  log_file=$Rf"/log"

  # Make the dot in file.
  # Note that the python code creates an empty f_string directory.
  # Also cd to where the ROMS executable lives.
  cd $LO_parent"/forcing/dot_in"
  source $HOME"/.bashrc"
  if [ $D = $D0 ] && [ $start_type = "new" ] ; then
    python ./make_dot_in.py $start_type $DD
    cd $R_parent/makefiles/cas1_baseM_tideramp
  else
    python ./make_dot_in.py continuation $DD
    cd $R_parent/makefiles/cas1_baseM
  fi

  # choose which cores to run on
  if [ $run_type = "forecast" ] ; then
    hf="hf72a"
  elif [ $run_type = "backfill" ] ; then
    hf="hf72b"
  fi

  # the actual ROMS run command
  if [ $HOME = "/Users/PM5" ] ; then # testing
    echo "/cm/shared/local/openmpi-ifort/bin/mpirun -np 72 -machinefile $hf oceanM $Rf/liveocean.in > $log_file &"
  elif [ $HOME == "/home/parker" ] ; then # the real thing
    /cm/shared/local/openmpi-ifort/bin/mpirun -np 72 -machinefile $hf oceanM $Rf/liveocean.in > $log_file &
    # Check that ROMS has finished successfully.
    checkfile=$log_file
    flag=0
    while [ $flag -ne 1 ]
    do
      sleep 30
      if [ -e $checkfile ]
      then
        last_line=$(tail -1 $checkfile)
        # example of good finish (note one space before ROMS)
        # ROMS/TOMS: DONE... Friday - November 21, 2014 -  5:54:03 PM
        # example of bad finish
        # ERROR: I/O related problem.
        if [[ $last_line = *"DONE"* ]]
        then
          flag=1
          echo $last_line
        elif [[ $last_line = *"ERROR"* ]]
        then
          echo "Problem with last line! "$checkfile
          echo $last_line
          flag=1 # exit this while loop
        else
          echo "assuming the run is still going"
          echo $last_line
        fi
      else
        # this executes if the checkfile does not exist
        echo "No checkfile yet"
      fi
    done
  fi
  
  # This function increments the day.
  # NOTE: it changes y, m, d, and D, even in the scope of this shell script!
  next_date $y $m $d
  
done

