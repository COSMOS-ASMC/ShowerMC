#!/bin/bash
 if [ $# != 2 ]; then
    echo "./procAll.sh inputfileDir outputFile"
    exit
 fi
#  rm -f $2
 for f in $1/*.hyb; do
    cosz=`awk '$1=="h" {print $9;exit}' $f `
    awk '$1=="t" {print $3, $8}' $f | ./tranFit 5  > temp
    xmax=`awk  '{print $1}' temp `
    axmax=`awk  '{print $2}' temp `
    
    awk '{print $1/cosz, $2/cosz, cosz}' cosz=$cosz  temp >> $2
 done
