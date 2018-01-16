#!/bin/bash

# test of the new get_env.sh code

# run the code to put the environment into a csv
if [ -e ../alpha/user_get_lo_info.sh ] ; then
  ../alpha/user_get_lo_info.sh
else
  ../alpha/get_lo_info.sh
fi

# and read the csv into active variables
while IFS=, read col1 col2
do
  # this sets the variable names and contents from the columns
  eval $col1="$col2"
  #
  # and this is just to have a look
  eval echo $col1" = "\$$col1
  #
done < ../alpha/lo_info.csv
