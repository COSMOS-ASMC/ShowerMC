#!/bin/bash

tmp=$0
SCRDIR=$(dirname $(readlink -f $tmp))
$SCRDIR/CleanUpTemporaryFiles.sh
if [ $? -ne 0 ]
then
    exit
fi

CURDIR=$(dirname $(readlink -f $SCRDIR))
if [ ! -d ${CURDIR}/build ]
then
    mkdir -p ${CURDIR}/build
fi

if [ $# -ge 1 ]
then
    BUILDDIR=$1
else
    BUILDDIR=${CURDIR}/build/ShowerMC
fi

if [ ! -d $BUILDDIR ]
then
    echo "argument error"
    echo "bash $0 LIBRARY_BUILD_DIRECTORY"
    echo "Do you want to set LIBRARY_BUILD_DIRECTORY as $BUILDDIR ? (y\n)"
    read yn
    if [ ${yn,,} == "y" ] || [ ${yn,,} == "yes" ]
    then
	mkdir -p $BUILDDIR
    else
	echo "abort."
	exit
    fi
fi

if [ ! $SHOWERMCTOP ] || [ ! $LIBLOFT ] || [ ! $COSMOSTOP ]
then
    source $SCRDIR/SetEnvironment.sh
fi

if [ ! $SHOWERMCTOP ] || [ ! $LIBLOFT ] || [ ! $COSMOSTOP ]
then
    echo Please set LIBLOFT and COSMOSTOP by yourself
    echo fatal error. abort.
    exit -1
fi

echo "SHOWERMCTOP is $SHOWERMCTOP"
echo "LIBLOFT     is $LIBLOFT"
echo "COSMOSTOP   is $COSMOSTOP"
echo "BUILDDIR    is $BUILDDIR"

### build library
cd $BUILDDIR

cmake ${SHOWERMCTOP}
if [ $? -ne 0 ];then
    echo "cmake error $?. abort."
    exit
fi
LANG=C make -j
if [ $? -ne 0 ];then
    LANG=C make
    if [ $? -ne 0 ];then
	echo "make error. abort."
	exit
    fi
fi

