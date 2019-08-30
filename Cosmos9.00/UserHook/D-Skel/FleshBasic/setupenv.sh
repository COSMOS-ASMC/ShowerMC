#!/bin/sh 
###  don't touch next  line:   test is needed for sge job.
test  $# -eq  0  &&    source ../confirm.inc

#   id to be used for the command file mentioned above. (must start with A-z)
EXECID=p15cos1.0
export EXECID
#------------may be you don't need touch next
source ../setupenv.sh
TOPDIR=$COSMOSTOP/UserHook/$ARENA
export TOPDIR
#   dir to store execution command (ssh)  or submit command (sge)
EXECDIR=$TOPDIR/FleshBasic/Exec
export EXECDIR
#     where param001 etc exist
PARAMDIR=$TOPDIR/FleshBasic/ParamDir
export PARAMDIR
#     where to put  files createdd by each host
#     /tmp/$USER/ is better
OUTDIR=$TOPDIR/Assemble/OutDir
#OUTDIR=/tmp/$USER
export OUTDIR
#     where to put error message
ERRDIR=$TOPDIR/FleshBasic/ErrDir
export ERRDIR

#
#  don't touch below.
#  if used from script, we skip the next lines.
if [ $# -eq  0 ] ; then
    confirm  $PARAMDIR
    confirm  $OUTDIR
    confirm  $ERRDIR
    confirm  $EXECDIR
    
    if [ -f Sparam ]; then
	primaryf=` awk '$1=="PRIMARYFILE" {i1=index($0, qm); i2=index(substr($0, i1+1), " "); print substr($0, i1+1,i2-1)}' qm="'" Sparam `
	echo "'"$primaryf"' seems to be the primary file used at SkelFlesh."
	echo "Now it is copied to this directory"
	cp ${COSMOSTOP}/UserHook/SkelFlesh/$primaryf ./
    else
	echo "There is no 'Sparam' file which should be a copy of "
	echo " parameter file used in the SkeleFlesh process"
	echo " **** make a copy of parameter file  as Sparam and change Job='newskel' "
        echo "into  Job='newflesh'"
	exit
    fi
fi
