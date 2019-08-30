#!/bin/sh
if [ $# -ne 5 ]; then
    echo "usage: conv2redT.sh MolierUnit Fai cosz input output"
    echo "convert (r,t) pairs into (r, reduced_time, dt, reduced_time')"
    echo "innput: (r,t) pair lines. r in MU t in ns."
    echo "output will have: r, redt, dt, redt'"
    exit 1
fi
awk -f $COSMSOSTOP/UserHook/mkLDD/Util/conv2redT.awk mu=$1 fai=$2 cosz=$3 $4 > $5
