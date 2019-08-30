#!/bin/sh
if [ $# != 1 ]; then
    echo "Arg.  1 is destination directory"
    exit
fi
#  set MCPU
source ../Smash/setupenv.sh $0
echo "MCPU=$MCPU"
#  we get OUTDIR,  GOUTDIR
source ../$FLESHDIR/setupenv.sh $0
GOUTDIR=${GOUTDIR:-$OUTDIR}
echo "OUTDIR=$OUTDIR  GOUTDIR=$GOUTDIR"
echo "EXECID=$EXECID"
source ./setupenvHyb.sh $0


nc=0

for f in $GOUTDIR/*.dat-r; do
   echo "procssing data $f"
# Usage  md2WebSysDat  input_dir file-body  output-dir 
  ./splitbyL $f $EXECID  $1/  
  nc=`expr $nc + 1`
  if [ $nc = $MCPU ]; then
      break
  fi
done
# restore stdin and close 3

if [ $MCPU -ne $nc ];  then
    echo "# of cpu's not enough; there might be running jobs"
    echo "$MCPU cpu's should exist but $nc is counted"
else
    echo "All events have been successfully assmebled to $1/"
fi
