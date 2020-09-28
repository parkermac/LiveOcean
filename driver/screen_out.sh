#!/bin/bash

# code to automate looking at screen output from a forecast forcing job

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

# defaults
gridname="cas6"
tag="v3"
ex_name="lo8b"
frc="ocn4"

while [ "$1" != "" ]; do
  case $1 in
    -g | --gridname )  shift
      gridname=$1
      ;;
    -t | --tag )  shift
      tag=$1
      ;;
    -f | --frc )  shift
      frc=$1
      ;;
  esac
  shift
done

gtag=$gridname"_"$tag

# do forecast only
ymd=$(date "+%Y.%m.%d")
f_string="f"$ymd
so=$LOo$gtag"/"$f_string"/"$frc"/Info/screen_output.txt"

cat $so