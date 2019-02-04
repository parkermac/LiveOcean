#!/bin/bash

# code to test tracker

ii=0

while [ $ii -le 30 ]
do
  echo $ii
  python pan_plot.py -lt merhab -pt P_tracks_ps_debug -test True
  ii=$[$ii + 1]
done