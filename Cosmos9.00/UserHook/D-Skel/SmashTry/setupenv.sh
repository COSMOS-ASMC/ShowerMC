#!/bin/sh 
## don't touch next line; test is needed for sge job  (may not be used from sge)
test  $# -eq  0  &&     source ../confirm.inc

#----------------
#      smash for FleshBasic or FleshHist
FLESHDIR=FleshHist
#FLESHDIR=FleshBasic
export FLESHDIR
#   number of cpu to be used in the distributed job : 
NCPU=900
export NCPU
#   MCPU <= NCPU;  if <, then fleshing is done only at MCPU  hosts and
#      cput power given in the hostlist must be equal.  
MCPU=45
export MCPU
#  For normal job with  NCPU=MCPU or NCPU > MCPU, this may be kept
#  as it is. It is computed as NCPU/MCPU inside.  However,
#  for the rescue job, you have to give NCPU=MCPU even if normal fleshing
#  job used NCPU > MCPU. In that case you have to give that ratio to this.
ENHANCE=20
export ENHANCE
MARGIN=4
export MARGIN
#    path to the binary skeleton file created by skeleton creation process
#    in SkelFlesh. You may use completely different path 
SKELETON=../$FLESHDIR/Skeleton
export SKELETON
#   directory to store smashed new  skeleton files to be created here
#   ( N smashed skeleton files ).  You may use any path
#   You must delete contents before smash proccess.

SKELDIR=../$FLESHDIR/SkelDir
export SKELDIR
#  base name of the smashed new skeleton files. 
#   If it is "Skeleton", then,  Skeleton001, Skeleton002 ... will be  put
#   in the  directory specified by SKELDIR

SKELNAME=Skeleton
export SKELNAME
#  If given,  it is a path to the file which contains a list of hosts
#  where distributed job will be performed.  The smash program
#	will read the relative cpu power of each host from that file ;
# 	partilces are	distributed so that the time for computation 
#	be almost equal on every host.
#	If not given, equal power is assumed.
HOSTLIST=../Host
export HOSTLIST
#      select Host from HOSTLIST randomy. if NCPU==MCPU, next is  not used.
THINHOSTLIST=../ThinHost
export THINHOSTLIST

#  don't touch below.
#  if used from script, we skip the next lines.
if [  $# -eq  0 ] ; then
    confirm  $SKELDIR
    dup=`sort -n $HOSTLIST | awk '{print $1}' | uniq -d`
    if [ "x$dup" != "x" ]; then
	echo "some numbers in $HOSTLIST duplicated"
	echo " they are:"
	echo "$dup"
	exit
    fi
fi
