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




#    00002  will be replaced by cpu# by execflesh_all.sh or _one.sh
NUMB=00002
export NUMB
source $COSMOSTOP/Scrpt/setarch.sh
#$  -e /home/Users/kasahara/Cosmos/UserHook/D-Skel/FleshHist/ErrDir/p20cos0.800E10-1-T1e-6x1e5-500keV-$JOB_ID-$HOSTNAME.00002.err
##  in the case of binary output, next can be omitted.
###  #$  -o /tmp/kasahara/p20cos0.800E10-1-T1e-6x1e5-500keV-$JOB_ID-$HOSTNAME.00002.dat
./DHflesh$ARCH < $PARAMDIR/param00002
