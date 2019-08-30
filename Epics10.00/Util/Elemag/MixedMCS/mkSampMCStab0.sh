#!/bin/bash
if [ $# != 2 ]; then
    cat << "EOF"
    Usage:  mkSampMCStab.sh  media-name  dir
    where
       media-name is such as Pb which is a medium name
          listed in $EPICSTOP/Data/BaseM/
       dir is the directory where the sampling table fileis stored
        The file name will be the same as media name
EOF
exit
fi
source $EPICSTOP/Scrpt/setarch.sh

echo $ARCH
#if [ ! -f mksmptab$ARCH.out ]; then
    make -f mkSmpTbl.mk
#fi
media=$1
dir=$2
echo $media | ./mksmptab$ARCH.out > $dir/$media
exit



