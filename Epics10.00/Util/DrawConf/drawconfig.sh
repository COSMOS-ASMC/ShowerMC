#!/bin/bash
work=/tmp/$USER/Work
if [ $# != 1 ]; then
    cat <<EOF
Usage: ./drawconfig.sh configfile-path
     where 
      configfile-path: path to a config
EOF
exit
fi
make clean 
make -f drawConfig.mk
mkdir -p $work
fullpath=$1
echo $1 > $work/.configrc
./drawconfig



