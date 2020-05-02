#!/bin/bash

# removes the extra hours (history files 0026---73) from all the
# output days of a user-specified run

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

# USE COMMAND LINE OPTIONS
#
# -g name of the grid [cascadia1, ...]
# -t name of the forcing tag [base, ...]
# -x name of the ROMS executable to use
# -0 start date: yyyymmdd
# -1 end date: yyyymmdd

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
  esac
  shift
done


gtag=$gridname"_"$tag
gtagex=$gtag"_"$ex_name

Rf=$roms"output/"$gtagex"/"$f_string"/"
