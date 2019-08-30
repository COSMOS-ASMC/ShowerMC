#!/bin/sh
#$  -S /bin/sh
dir_bin=/home/Users/kasahara/Cosmos/UserHook/FirstKiss
dir_out=/home/Users/kasahara/Cosmos/UserHook/FirstKiss
COSMOSTOP=/home/Users/kasahara/Cosmos
export COSMOSTOP
export dir_bin
export dir_out
source /home/Users/sgeadmin/default/common/settings.sh
source  /Loft1/Intel/ifc/bin/ifortvars.sh
#$ -o /home/Users/kasahara/Cosmos/UserHook/FirstKiss/myout
#$ -e /home/Users/kasahara/Cosmos/UserHook/FirstKiss/myerr
pwd
cd $dir_bin
pwd
echo ls -trl
echo
echo  start cosmos
echo ===============
./cosmosPCLinuxIFC < param
date
pwd
ls -trl
date
date
