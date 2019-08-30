#!/bin/sh
#PBS -q B
#      next;  stdout and stderr --> same file
# #PBS -j oe
#PBS -e /tibet/work/kasahara/ErrDir2/${PBS_JOBID}.err

COSMOSTOP=/tibet/work/kasahara/Cosmos7.52
export COSMOSTOP
cd $COSMOSTOP/UserHook/mkLDD

source  /tibet/work/kasahara/intel/ifc/ifortvars.sh intel64
#source  /opt/intel/cc/9.1.052/bin/iccvars.sh
source  /tibet/work/kasahara/intel/icc/iccvars.sh intel64
#source  /tibet/work/kasahara/intel/ifc/bin/ifortvars.sh 
if [ ! "${HOST}" ];then
    export HOST=$HOSTNAME
fi

#eeee
#./cosmosPCLinuxIFC < param  >/dev/null 
./cosmosPCLinuxIFC64 < param 
##  sync -avz  /tmp/$USER/EXECID-*.__.* /tibet/work/$USER/$ARENA/

