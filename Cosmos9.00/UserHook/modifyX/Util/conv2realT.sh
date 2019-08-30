#!/bin/sh
if [ $# -ne 5 ]; then
    echo "usage: conv2realT.sh MolierUnit Fai cosz input output"
    echo "output will have: r, t, dt, t+dt0"
    echo "If last value <0, indicate numerical error for T~0"
    exit 1
fi
awk -f $COSMOSTOP/UserHook/mkLDD/Util/conv2realT.awk mu=$1 fai=$2 cosz=$3 $4 > $5


