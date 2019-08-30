#!/bin/csh

if( $#argv != 1) then
    echo Usage: ./sortptcls.csh  timesorted_data
    exit
endif

awk -f sortptcls.awk 
awk  '$4==3 {print $1, $2, $3}' $1 > muon.dat
awk '$4 > 32 {print $1, $2, $3}' $1 > other.dat

