#!/bin/bash
if [ ! -s Work/tList ] ; then
    echo "need"
else
    echo "other"
fi

echo $#
nc=1
let nc=nc+1
echo $nc
let nc=$nc+1
echo $nc
let nc++
echo $nc
ARCH=`./setarch.sh`
echo $ARCH
