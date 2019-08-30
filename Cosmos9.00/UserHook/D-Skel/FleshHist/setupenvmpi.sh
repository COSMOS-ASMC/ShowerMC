#!/bin/bash 
###  don't touch next  line.   test is needed for sge job.
test  $# -eq  0  &&    source ../confirm.inc
n=0
(cd ./SmSkelDir; 
  rm -f rank-*;
  for f in *; do
    echo "smashed skeleton $f is being linked to rank-$n"
    ln -s $f rank-$n
    n=`expr $n + 1`
  done
)
if [ $# -eq  0 ] ; then
    confirm  /bgwork0/sctasim/kasahara/FDD2
fi
echo "In Sparam, UserHookc must be modified so that the first value"
echo "is /bghome/sctasim/kasahara/Cosmos/UserHook/D-Skel/FleshHist/SmSkelDir/rank-#"
sel=1
while [ $sel -ne 0 ]; do
    echo " "
    echo "select number"
    echo "0- quit"
    echo "1- edit Sparam"
    echo "2- edit jsfFlesh"
    echo "3- make executable (make clean;make)"
    echo "4- see job class status(bllclass)"
    echo "5- submit job by bllsubmit jsfFlesh"
    read sel

    if [ $sel -eq 1 ]; then
	emacs ./Sparam
    elif [ $sel -eq 2 ]; then
	emacs ./jsfFlesh
    elif [ $sel -eq 3 ]; then
	make clean;make
    elif [ $sel -eq 4 ]; then
	bllclass
    elif [ $sel -eq 5 ]; then
	bllsubmit jsfFlesh
	echo "status: see ErrDir/err"
	sleep 3
    elif [ $sel -eq 0 ]; then
	break;
    fi
done
