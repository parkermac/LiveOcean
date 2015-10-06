#!/bin/bash

# This fills in gaps in our WRF inventory.
# Designed to be run on fjord

dir0="/pmraid3/darr/tstwrf/tmpwrf/"

d1="2014010612"
dmissing="2014010700"
d2="2014010712"

cd $dir0

echo "d1"
ls $d1
echo "d2"
ls $d2
