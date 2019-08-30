#!/bin/sh
#   sge
#$ -S /bin/sh 
ARENA=D-Skel
COSMOSTOP=~/Cosmos
export COSMOSTOP
cd $COSMOSTOP/UserHook/$ARENA/FleshHist
source ./setupenv.sh $0
source ../Smash/setupenv.sh $0
#bbbbb  You may have to change these two lines depending
######## on your env.
source  /home/Users/sgeadmin/SGE/default/common/settings.sh 
source  /Loft1/Intel/ifc/bin/ifortvars.sh
#eeee




#    __  will be replaced by cpu# by execflesh_all.sh or _one.sh
NUMB=__
export NUMB
source $COSMOSTOP/Scrpt/setarch.sh
#$  -e ERRDIR/EXECID-$JOB_ID-$HOSTNAME.__.err
##  in the case of binary output, next can be omitted.
###  #$  -o OUTDIR/EXECID-$JOB_ID-$HOSTNAME.__.dat
./DHflesh$ARCH < $PARAMDIR/param__
