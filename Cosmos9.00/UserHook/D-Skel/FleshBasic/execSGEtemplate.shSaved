#!/bin/sh
#   sge
#$ -S /bin/sh 

COSMOSTOP=~/Cosmos
export COSMOSTOP
source ../setupenv.sh 
source $COSMOSTOP/UserHook/$ARENA/FleshBasic/setupenv.sh $0

#bbbbb  You may have to change these two lines depending
######## on your env.
source  /home/Users/sgeadmin/default/common/settings.sh
source  /Loft1/Intel/ifc/bin/ifortvars.sh
#eeee


cd $TOPDIR/FleshBasic
source ./setupenv.sh $0
source ../Smash/setupenv.sh $0

#    __  will be replaced by cpu# by execflesh_all.sh or _one.sh
NUMB=__
export NUMB
source $COSMOSTOP/Scrpt/setarch.sh
#$  -e ERRDIR/EXECID-$JOB_ID-$HOSTNAME.__.err
#$  -o OUTDIR/EXECID-$JOB_ID-$HOSTNAME.__.dat
./flesh$ARCH < $PARAMDIR/param__


