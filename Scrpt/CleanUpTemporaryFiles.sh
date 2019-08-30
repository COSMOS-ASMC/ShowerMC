#!/bin/bash

tf=`find ./ -name "cmake_install.cmake" -o -name "temp*.[fF]*" -o -name "*.mod" -o  -name "*.o" -o -name "#*" -o -name ".#*" -o -name "*~" | grep -v lib/`
if [ "$tf" ]
then
    echo Remove following files before compiling
    echo $tf
    echo "Do you want to remove them now ? (y\n)"
    read yn
    if [ $yn == "y" ] || [ $yn == "yes" ] || [ $yn == "YES" ] || [ $yn == "Yes" ]
    then
	rm -rf $tf
    else
	echo "abort."
	exit -1
    fi
fi
