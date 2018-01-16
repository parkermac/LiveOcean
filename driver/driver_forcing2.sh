#!/bin/bash

# This runs the code to create forcing for one or more days,
# for any of the types of forcing, allowing for either a
# forecast or backfill over a range.

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

. $LO"driver/common.lib"

# USE COMMAND LINE OPTIONS
#
# -g name of the grid [cascadia1]
# -t name of the forcing tag [base]
# -x name of the ROMS executable to use; only needed for post processing [lobio5]
# -f forcing type [atm, ocn, riv, tide, azu, low_pass, etc.]
# -r run type [forecast, backfill]
#  if backfill, then you must provide two more arguments
# -0 start date: yyyymmdd
# -1 end date: yyyymmdd
# -nc do not remake forcing if it already exists [no argument]
#
# example call to do backfill:
# ./driver_forcing1.sh -g cascadia1 -t base -f atm -r backfill -0 20140201 -1 20140203
#
# example call to do forecast:
# ./driver_forcing1.sh -g cascadia1 -t base -f atm -r forecast
#
# example call push to azure:
# ./driver_forcing1.sh -g cascadia1 -t base -x lobio5 -f azu -r backfill -0 20140201 -1 20140203

ex_name="placeholder"
clobber_flag=1 # the default (1) is to clobber, unless the -nc argument is used
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
    -nc | --do_not_clobber )
      clobber_flag=0
      ;;
  esac
  shift
done

if [ $run_type = "forecast" ] ; then
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

keep_going=1 # 1 => keep going, 0 => stop the driver
while [ $D -le $D1 ] && [ $keep_going -eq 1 ]
do
  # manipulate the string D to insert dots, using the syntax:
  #    substring = ${string:start_index:count}
  # and the index starts from 0
  DD=${D:0:4}.${D:4:2}.${D:6:2}
  
  f_string="f"$DD
  LOog=$LOo$gtag"/"
  LOogf=$LOog$f_string"/"
  LOogf_f=$LOogf$frc"/"
  echo $LOogf_f
  echo $(date)
  LOogf_fi=$LOogf_f"Info/"
  LOogf_fd=$LOogf_f"Data/"
  
  # this makes all parent directories if needed
  mkdir -p $LOogf
  
  checkfile=$LOogf_fi"process_status.csv"
  
  already_done_flag=0
  # check to see if the job has already completed successfully
  if [ -f $checkfile ] && [ $clobber_flag -eq 0 ]; then
    if grep -q "result,success" $checkfile; then
      echo "- No action needed: job completed successfully already."
      already_done_flag=1
    else
      echo "- Trying job again."
    fi
  fi
  
  if [ $already_done_flag -eq 0 ] || [ $clobber_flag -eq 1 ]; then
    
    if [ -d $LOogf_f ]
    then
      rm -rf $LOogf_f
    fi
    mkdir $LOogf_f
    mkdir $LOogf_fi
    mkdir $LOogf_fd
    
    # Make the forcing.
    cd $LO"forcing/"$frc
    if [ -e $HOME"/.bashrc" ] ; then
      source $HOME"/.bashrc"
    fi
    if [ -e $HOME"/.bash_profile" ] ; then
      source $HOME"/.bash_profile"
    fi
    if [ -e $HOME"/.profile" ] ; then
      source $HOME"/.profile"
    fi
    
    
    python ./make_forcing_main.py -g $gridname -t $tag -f $frc -r $run_type -d $DD -x $ex_name > $LOogf_fi"screen_out.txt" &

    # Check that the job has finished successfully.
    PID1=$!
    wait $PID1
    echo "job completed for" $f_string
    echo $(date)
  
    # check the checkfile to see if we should continue
    if grep -q "result,success" $checkfile ; then
      echo "- Job completed successfully."
      # will continue because we don't change keep_going
    else
      echo "- Something else happened."
      keep_going=0
      # stop the driver
    fi
  
  fi # end of already_done_flag test
  
  # This function increments the day.
  # NOTE: it changes y, m, d, and D, even in the scope of this shell script!
  next_date $y $m $d
  
done
