#!/bin/sh
COSMOSTOP=~/Cosmos
export COSMOSTOP
# source  /Loft1/Intel/ifc/bin/ifortvars.sh
cd $COSMOSTOP/UserHook/mkLDD
( ./cosmosPCLinuxIFC < param >  /dev/null ) >& ./ErrDir/$HOSTNAME.$$.err &













