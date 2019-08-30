#!/bin/bash

if [ $# -lt 1 ] || [ ! -d $1 ]
then
    echo "argument error"
    echo "bash $0 TARGET_DIRECTORY [BUILD_DIRECTORY]"
    echo abort.
    exit
fi

if [ ! $SHOWERMCTOP ] || [ ! $LIBLOFT ] || [ ! $COSMOSTOP ]
then
    source $(dirname $0)/SetEnvironment.sh
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

tmp=$1
TARGETDIR=$(readlink -f $tmp)
TARGETNAME=$(basename $TARGETDIR)
if [ $# -eq 1 ];then
    dir=${SHOWERMCTOP}/build/${TARGETNAME}
    if [ ! -d $dir ]
    then
	echo "Do you want to set BUILD_DIRECTORY as $dir ? (y\n)"
	read yn
	if [ $yn == "y" ] || [ $yn == "yes" ] || [ $yn == "YES" ] || [ $yn == "Yes" ]
	then
	    mkdir -p $dir
	else
	    echo "abort."
	    exit
	fi
    fi
else
    dir=$2
fi


### build library
cd $dir

cmake ${TARGETDIR}
if [ $? -ne 0 ];then
    echo "cmake error. abort."
fi
LANG=C make -j
if [ $? -ne 0 ];then
    LANG=C make
    if [ $? -ne 0 ];then
	echo "make error $?. abort."
	exit
    fi
fi

