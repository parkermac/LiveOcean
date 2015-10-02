#!/bin/bash

# test of input parameter code
# 9/16/2015 PM

# RESULT: this works correctly, parsing the input regardless of order

while [ "$1" != "" ]; do
  case $1 in
    -a )  shift
      a=$1
      ;;
    -b )  shift
      b=$1
      ;;
    -c )  shift
      c=$1
  esac
  shift
done

echo "a "$a
echo "b "$b
echo "c "$c