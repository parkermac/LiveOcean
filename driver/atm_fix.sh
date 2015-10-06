#!/bin/bash

# This fills in gaps in our WRF inventory.
# Designed to be run on fjord

dir00="/pmraid3/darr/tstwrf/tmpwrf/"
dmissing="TEST_2014010700"
mkdir $dmissing
cd $dir00

dir1="2014010612" # get f12 to f24 and rename as f00 to f12 in dmissing
nout=0
while [ $nout -le 12 ]
do
  
  if [ $nout -le 9 ] ; then
    pad="0"
  else
    pad=""
  fi
  outfile=$dmissing"/wrfout.ocean_d2."$dmissing".f"$pad$nout".0000"
  
  nin=$[$nout + 12]
  if [ $nin -le 9 ] ; then
    pad="0"
  else
    pad=""
  fi
  infile=$dir1"/wrfout.ocean_d2."$dir1".f"$pad$nin".0000"
  
  echo "copying "$infile" to "$outfile
  
  cp $infile $outfile
  
  nout=$[$nout + 1]
  
done

dir1="2014010712" # get f01 to f12 and rename as f13 to f24 in dmissing
nout=13
while [ $nout -le 24 ]
do
  
  if [ $nout -le 9 ] ; then
    pad="0"
  else
    pad=""
  fi
  outfile=$dmissing"/wrfout.ocean_d2."$dmissing".f"$pad$nout".0000"
  
  nin=$[$nout - 12]
  if [ $nin -le 9 ] ; then
    pad="0"
  else
    pad=""
  fi
  infile=$dir1"/wrfout.ocean_d2."$dir1".f"$pad$nin".0000"
  
  echo "copying "$infile" to "$outfile
  
  #cp $infile $outfile
  
  nout=$[$nout + 1]
  
done
#wrfout.ocean_d2.2014043012.f00.0000
