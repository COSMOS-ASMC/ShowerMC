#!/bin/sh
#PBS -q A
#      next;  stdout and stderr --> same file
## #PBS -j oe
#PBS -o /tibet/work/kasahara/EpicsData/Invade/Genpi0/qgsICR-${PBS_JOBID}.dat
#PBS -e /tibet/work/kasahara/ErrDir/${PBS_JOBID}.err

COSMOSTOP=/tibet/work/kasahara/Develop/Cosmos7.53
export COSMOSTOP
EPICSTOP=/tibet/work/kasahara/Develop/Epics9.00
export EPICSTOP

source  /tibet/work/kasahara/intel/ifc/ifortvars.sh intel64
#source  /opt/intel/cc/9.1.052/bin/iccvars.sh
source  /tibet/work/kasahara/intel/icc/iccvars.sh intel64
#source  /tibet/work/kasahara/intel/ifc/bin/ifortvars.sh 
if [ ! "${HOST}" ];then
    export HOST=$HOSTNAME
fi

cd $EPICSTOP/Util/Gencol
#eeee
#./cosmosPCLinuxIFC < param  >/dev/null 
./GencolPCLinuxIFC64 < param | awk '$1==4 && $3==0 && $4>600' 
##  sync -avz  /tmp/$USER/EXECID-*.__.* /tibet/work/$USER/$ARENA/
