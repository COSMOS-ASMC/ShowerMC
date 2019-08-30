#!/bin/bash

dir=$(readlink -f $(dirname $BASH_SOURCE)/..)
    export SHOWERMCTOP=$dir
if [ $# -eq 2 ] && [ -d $1 ] && [ -d $2 ]
then
    export LIBLOFT=`readlink -f $1`
    export COSMOSTOP=`readlink -f $2`
    return
fi

echo If you want to set environment variable manually, 
echo \# source $BASH_SOURCE LIBLOFT_TOP_DIRECTORY COSMOS_TOP_DIRECTORY

TMP=$(ls -dt $(find $dir -type d -iname "libloft*"))
for i in $TMP
do
    if [ -d $i/Header ]
    then
	j=$(readlink -f $i)
	echo "Do you want to set LIBLOFT as $j ? (y\n)"
	read yn
	if [ ${yn,,} == "y" ] || [ ${yn,,} == "yes" ]
	then
	    export LIBLOFT=$j
	    break
	fi
    fi
done

TMP=$(ls -dt $(find $dir -type d -iname "cosmos*"))
for i in $TMP
do
    if [ -d $i/cosmos ]
    then
	j=$(readlink -f $i)
	echo "Do you want to set COSMOSTOP as $j ? (y\n)"
	read yn
	if [ ${yn,,} == "y" ] || [ ${yn,,} == "yes" ]
	then
	    export COSMOSTOP=$j
	    break
	fi
    fi
done

if [ ! $LIBLOFT ] || [ ! $COSMOSTOP ]
then
    echo Please set LIBLOFT and COSMOSTOP by yourself
    echo fatal error. abort.
    return
fi

export PATH=${LIBLOFT}/Scrpt:${COSMOSTOP}/Scrpt:$PATH
