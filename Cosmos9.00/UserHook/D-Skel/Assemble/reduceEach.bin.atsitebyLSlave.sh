#!/bin/sh
#   this is used by reduceEach.bin.atsite.csh
#  .nrfai for all events
ARENA=D-Skel
cd  $COSMOSTOP/UserHook/$ARENA/Assemble
NRFAIFILE=$3
export NRFAIFILE

#  .dat files
datfiles=/tmp/${USER}/*$1.$2.dat

rm -f /tmp/${USER}/*$1.$2.dat-r
rm -f  /tmp/${USER}/*$1.$2.nrfai-r
#-----------------
source ../FleshHist/setupenv.sh  $0
for  f in $datfiles
do
     echo $f >> error
     DATFILE=$f
     export DATFILE
     outfile=$f"-r"
     if [  -f $f ]; then
         $COSMOSTOP/UserHook/$ARENA/Assemble/reduceEachSize.binbyL.PCLinuxIFC > $outfile 2>> error  &
#           mv -f $outfile $f-r
    fi
done

