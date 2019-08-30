#!/bin/sh
#   sge
#$ -S /bin/sh 

COSMOSTOP=~/Cosmos
export COSMOSTOP
source $COSMOSTOP/UserHook/DisParaForTA/FleshHist/setupenv.sh $0
source $COSMOSTOP/UserHook/DisParaForTA/Smash/setupenv.sh $0

#bbbbb  You may have to change these two lines depending
######## on your env.
source  /home/Users/sgeadmin/default/common/settings.sh 
source  /Loft1/Intel/ifc/bin/ifortvars.sh
#eeee


cd $TOPDIR/FleshHist

#    __  will be replaced by cpu# by execflesh_all.sh or _one.sh
NUMB=__
export NUMB
source $COSMOSTOP/Scrpt/setarch.sh
#$  -e ERRDIR/EXECID-$JOB_ID-$HOSTNAME.__.err
##  in the case of binary output, next can be omitted.
###  #$  -o OUTDIR/EXECID-$JOB_ID-$HOSTNAME.__.dat
x=` echo $OUTDIR | grep "/tmp/" `
if [ x$x = "x" ]; then
  ./DHflesh$ARCH < $PARAMDIR/param__
else
#  ./DHflesh$ARCH < $PARAMDIR/param__  ; rsync -e ssh -avz $OUTDIR/ tasim502:$OUTDIR/
  ./DHflesh$ARCH < $PARAMDIR/param__ 
fi

