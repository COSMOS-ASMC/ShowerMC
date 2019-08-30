#!/bin/bash

CURDIR=`pwd`

SCRDIR=$(dirname $(readlink -f $0))
$SCRDIR/CleanUpTemporaryFiles.sh
if [ $? -ne 0 ]
then
    exit
fi

if [ "$1" = "keep" ]
then
    log=$SCRDIR/../tmp/makelog.txt
else
    log=/dev/null
fi

if [ ! $LIBLOFT ] || [ ! $COSMOSTOP ]
then
    source $SCRDIR/SetEnvironment.sh
fi

if [ ! $LIBLOFT ] || [ ! $COSMOSTOP ]
then
    echo Please set LIBLOFT and COSMOSTOP by yourself
    echo fatal error. abort.
    exit -1
fi

echo LIBLOFT   is $LIBLOFT
echo COSMOSTOP is $COSMOSTOP

deadlink=`find $SCRDIR/.. -xtype l`
if [ "$deadlink" ]
then
    echo WARNING : deadlink exist
    ls -l $deadlink
    echo
fi

cd $LIBLOFT
LANG=C make | tee $log
cd $COSMOSTOP
LANG=C make | tee -a $log

