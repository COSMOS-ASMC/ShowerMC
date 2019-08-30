#!/bin/bash
 if [ $# != 2 ]; then
    echo "./getmax.sh inputfileDir outputFile"
    exit
 fi
 rm -f $2
 for f in $1/*.hyb; do
    cosz=`awk '$1=="h" {print $9}' $f `
    awk '$1=="t" {print $3, $8}' $f | ./tranFit 3 | awk '{print $1/cosz, $2/cosz, cosz}' cosz=$cosz  >> $2
 done