#!/bin/sh -f
# ssh
source ../setupenv.sh
cd $COSMOSTOP/UserHook/$ARENA/FleshBasic/
source ./setupenv.sh $0
source ../Smash/setupenv.sh $0
#    __ will be replaced by cpu# by execflesh_all or _one
NUMB=__
export NUMB
source $COSMOSTOP/Scrpt/setarch.sh
(./flesh$ARCH < $PARAMDIR/param__  > $OUTDIR/$EXECID-$HOST.__.dat ) >& $ERRDIR/$EXECID-$HOST.__.err & 
