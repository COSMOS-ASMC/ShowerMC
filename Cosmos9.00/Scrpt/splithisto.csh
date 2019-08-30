#!/bin/csh -f
if( $#argv != 2 )  then
    echo "Usage: "$0 " maindir ascii-histfile"
    exit
endif
gawk -f $COSMOSTOP/Scrpt/splithisto.awk maindir=$1 $2




