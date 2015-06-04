#!/bin/bash

# This runs the code to create forcing for one or more days,
# for any of the types of forcing,
# allowing for either a forecast or backfill.

# NOTE: must be run from fjord.

# set a path and connect to a library of functions
if [ $HOME = "/Users/PM5" ] ; then
  LO_parent="/Users/PM5/Documents/LiveOcean"
elif [ $HOME = "/home/parker" ] ; then
  LO_parent="/data1/parker/LiveOcean"
fi
. $LO_parent"/driver/common.lib"

# COMMAND LINE INPUTS
#
# 1. which force (atm, ocn, riv, tide)
# 2. forecast or backfill
#  if backfill, then you must provide two more arguments
#  3. start date: yyyymmdd
#  4. end date: yyyymmdd
#
# example call to do backfill:
# ./driver_forcing.sh atm backfill 20140201 20140203 > log_forcing_atm &
# example call from cron:
# ./driver_forcing.sh atm forecast

frc=$1
run_type=$2

# check command line inputs
if [ $# -eq 2 ] ; then
  if [ $frc = "atm" ] || [ $frc = "ocn" ] || [ $frc = "riv" ] || [ $frc = "tide" ] || [ $frc = "azu" ] || [ $frc = "bio" ]
  then
    # do forecast
    ymd0=$(date "+%Y%m%d")
    ymd1=$ymd0
  else
    echo 'Bad input arguments (two)'
    exit
  fi
elif [ $# -eq 4 ] ; then
  if [ $frc = "atm" ] || [ $frc = "ocn" ] || [ $frc = "riv" ] || [ $frc = "tide" ] || [ $frc = "azu" ] || [ $frc = "bio" ] && [ $run_type = "backfill" ]
  then
    # do backfill
    ymd0=$3 # first date
    ymd1=$4 # last date
  else
    echo 'Bad input arguments (four)'
    exit
  fi
else
  echo 'Need 2 or 4 input arguments'
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
  echo $DD
  
  f_string="f"$DD
  LOo=$LO_parent"_output"
  LOog=$LOo"/"$gtag
  LOogf=$LOog"/"$f_string
  echo $LOogf
  LOogf_f=$LOogf"/"$frc
  LOogf_fi=$LOogf_f"/Info"
  LOogf_fd=$LOogf_f"/Data"
  # tedious... unnecessary?
  mkdir $LOo
  mkdir $LOog
  mkdir $LOogf
  if [ -d $LOogf_f ]
  then
    rm -rf $LOogf_f
  fi
  mkdir $LOogf_f
  mkdir $LOogf_fi
  mkdir $LOogf_fd
  
  # Make the forcing.
  cd $LO_parent"/forcing/"$frc
  source $HOME"/.bashrc"
  if [ $frc == "azu" ] ; then
    /home/parker/anaconda/bin/python ./make_forcing_main.py $frc $run_type $DD > $LOogf_fi"/screen_out.txt" &
  else
    python ./make_forcing_main.py $frc $run_type $DD > $LOogf_fi"/screen_out.txt" &
  fi

  # wait a bit to allow main to get rid of the output
  sleep 5
   
  # Check that the worker has finished successfully.
  # If it has then this file will exist, and its last line should read:
  # result,success
  checkfile=$LOogf_fi"/process_status.csv"
  
  flag=0
  while [ $flag -ne 1 ]
  do
    sleep 5
    if [ -e $checkfile ] ; then
      last_line=$(tail -1 $checkfile)
      if [ "$last_line" = "result,success" ] ; then
        flag=1
      else
        echo "Problem with last line! "$checkfile
        flag=1 # exit this while loop anyway
      fi
    fi
  done
  
  # print the contents of checkfile
  echo " "
  while read line
  do
    echo $line
  done <$checkfile
  echo " "
  
  # This function increments the day.
  # NOTE: it changes y, m, d, and D, even in the scope of this shell script!
  next_date $y $m $d
  
done
