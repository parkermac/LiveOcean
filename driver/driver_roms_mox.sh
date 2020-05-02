#!/bin/bash

# This runs ROMS for one or more days, allowing for either
# a forecast or backfill.
#
# Designed to be run from MOX
# and depends on other drivers having been run first on BOILER
#
# Designed to allow more flexible command line specification, and
# less information buried in the makefile programs.  This is about the
# -np and -N flags which specify the total number of cores, and the
# cores per node.  For the current mox environment, acceptable choices are
# -np 64 -N 32
#   or
# -np 196 -N 28
# i.e. np has to be an even multiple of N, and N has to be <= the number
# of nodes of that size that I own.

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
# -g name of the grid [cas5]
# -t name of the forcing tag [v3]
# -x name of the ROMS executable to use [lo8]
# -np number of cores to use [196 for N=28, 64 for N=32]
# -N cores per node [28 or 32 on mox]
# -s new or continuation
# -r forecast or backfill
#  if backfill, then you must provide two more arguments:
# -0 start date: yyyymmdd
# -1 end date: yyyymmdd
#
# example calls:
# ./driver_roms_mox.sh -g cas4 -t v2 -x lo6biom -np 64 -N 32 -s continuation -r forecast
# ./driver_roms_mox.sh -g cas5 -t v3 -x lo8 -np 196 -N 28 -s new -r backfill -0 20170101 -1 20170101
# ./driver_roms_mox.sh -g cas5 -t v3 -x lo8 -np 196 -N 28 -s continuation -r backfill -0 20170102 -1 20170131

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
    -np | --np_num )  shift
      np_num=$1
      ;;
    -N | --cores_per_node )  shift
      cores_per_node=$1
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
  echo "********** ROMS Driver *******************"
  echo "  blow ups = " $blow_ups
  
  # manipulate the string D to insert dots, using the syntax:
  # substring = ${string:start_index:count} and the index starts from 0
  DD=${D:0:4}.${D:4:2}.${D:6:2}

  f_string="f"$DD
  Rf=$roms"output/"$gtagex"/"$f_string"/"
  log_file=$Rf"log.txt"
  
  # get the f_string for two days before, so those can be deleted
  two_days_before $y $m $d
  DDP=${DPP:0:4}.${DPP:4:2}.${DPP:6:2}
  f_string_previous="f"$DDP
  # echo "TESTING: f_string = "$f_string", f_string_previous="$f_string_previous

  # Run make_dot_in.py, which creates an empty f_string directory,
  # and then cd to where the ROMS executable lives.

  cd $LO"forcing/dot_in/"$gtagex
  
  if [ -e $HOME"/.bashrc" ] ; then
    source $HOME"/.bashrc"
  elif [ -e $HOME"/.bash_profile" ] ; then
    source $HOME"/.bash_profile"
  elif [ -e $HOME"/.profile" ] ; then
    source $HOME"/.profile"
  fi

  if [ $D == $D0 ] && [ $start_type == "new" ] ; then
    python ./make_dot_in.py -g $gridname -t $tag -s new -r $run_type -d $DD -x $ex_name -np $np_num -bu $blow_ups
    sleep 30
    cd $roms"makefiles/"$ex_name
  else
    python ./make_dot_in.py -g $gridname -t $tag -s continuation -r $run_type -d $DD -x $ex_name -np $np_num -bu $blow_ups
    sleep 30
    cd $roms"makefiles/"$ex_name
  fi

  # run ROMS
  if [ $lo_env == "pm_mac" ] ; then # testing
    python make_back_batch.py -xp $Rf -np $np_num -N $cores_per_node -x $ex_name
    echo " -- just testing --"
    keep_going=1
    
  elif [ $lo_env == "pm_mox" ] ; then
    
    remote_dir="parker@boiler.ocean.washington.edu:/data1/parker/"
    local_dir="/gscratch/macc/parker/"
    
    # 1. Get today's forcing files from boiler (just once per day):
    if [ $blow_ups -eq 0 ] ; then
      scp -r $remote_dir"LiveOcean_output/"$gtag/$f_string $local_dir"LiveOcean_output/"$gtag
      PID1=$!
      wait $PID1
      echo "Forcing files transferred for "$f_string
    fi
    
    # 2. Run ROMS
    
    python make_back_batch.py -xp $Rf -np $np_num -N $cores_per_node -x $ex_name

    sbatch -p macc -A macc lo_back_batch.sh &
    
    # check the log_file to see if we should continue
    keep_checking_log=1
    while [ $keep_checking_log -eq 1 ]
    do
      sleep 30
      if [ -e $log_file ] ; then
        if grep -q "Blowing-up\|BLOWUP" $log_file ; then
          echo "- Run blew up!"
          blow_ups=$(( $blow_ups + 1 )) #increment the blow ups
          keep_checking_log=0
          if [ $blow_ups -le 7 ] ; then
            keep_going=1
          else
            keep_going=0
          fi
        elif grep -q "ERROR" $log_file ; then
          echo "- Run had an error."
          keep_going=0
          keep_checking_log=0
        elif grep -q "ROMS/TOMS: DONE" $log_file ; then
          echo "- Run completed successfully."
          keep_going=1
          blow_ups=0
          keep_checking_log=0
          
          # 3. Copy ROMS files to boiler:
          scp -r $local_dir"LiveOcean_roms/output/"$gtagex/$f_string $remote_dir"LiveOcean_roms/output/"$gtagex 
          PID1=$!
          wait $PID1
          echo "ROMS output files transferred for "$f_string
          
          # 4. Delete forcing and ROMS files from mox from two days before.
          echo "Deleting history and forcing files for "$f_string_previous
          rm -rf $local_dir"LiveOcean_output/"$gtag/$f_string_previous
          rm -rf $local_dir"LiveOcean_roms/output/"$gtagex/$f_string_previous
                    
        fi
      fi
    done
    echo "run completed for "$f_string
    
  fi
  
  echo $(date)
  echo " "

  if [ $blow_ups -eq 0 ] ; then
    # This function increments the day.
    # NOTE: it changes y, m, d, and D, even in the scope of this shell script!
    next_date $y $m $d
    blow_ups=0 # reset the blow up counter
  fi

done # end of while loop