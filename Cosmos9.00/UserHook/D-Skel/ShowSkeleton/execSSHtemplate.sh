#!/bin/sh  -f
# ssh
cd $COSMOSTOP/UserHook/DisParaForTA2/ShowSkeleton
source ./setupenv.sh $0
source ../Smash/setupenv.sh $0
#   __  will be replaced by cpu# by execflesh_all or _one
NUMB=__
export NUMB
source $COSMOSTOP/Scrpt/setarch.sh
#  when ascii output
(./DHflesh$ARCH < $PARAMDIR/param__  > $OUTDIR/$EXECID-$HOST.__.dat ) >& $ERRDIR/$EXECID-$HOST.__.err & 
#  when bin output
# ./DHflesh$ARCH < $PARAMDIR/param__   >& $ERRDIR/$EXECID-$HOST.__.err & 
