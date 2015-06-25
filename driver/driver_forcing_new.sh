#!/bin/bash

# This runs the code to create forcing for one or more days,
# for any of the types of forcing,
# allowing for either a forecast or backfill.

# This NEW version (6/24/2015) is designed to allow running a different ROMS
# executable using a command line argument.

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
# -x name of the ROMS executable to use (only needed if forcing is "azu") [lo1, ...]
# -f forcing type [atm, ocn, riv, tide, azu]
# -r run type [forecast, backfill]
#  if backfill, then you must provide two more arguments
# -0 start date: yyyymmdd
# -1 end date: yyyymmdd
#
# example call to do backfill:
# ./driver_forcing_new.sh -f atm -r backfill -0 20140201 -1 20140203
#
# example call to do forecast:
# ./driver_forcing_new.sh -f atm -r forecast
#
# example call push to azure:
# ./driver_forcing_new.sh -x lo1 -f azu -r backfill -0 20140201 -1 20140203
#
# you can also use long names like --ex_name instead of -x
#
# note that the --ex_name (or -x) parameter is only needed when calling azu.

while [ "$1" != "" ]; do
  case $1 in
    -x | --ex_name )  shift
      ex_name=$1
      ;;
    -f | --frc )  shift
      frc=$1
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
    /home/parker/anaconda/bin/python ./make_forcing_main.py $frc $run_type $DD -x $ex_name > $LOogf_fi"/screen_out.txt" &
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
