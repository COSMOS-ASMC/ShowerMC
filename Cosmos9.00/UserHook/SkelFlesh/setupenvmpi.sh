#!/bin/sh 
## don't touch next line; test is needed for sge job  (may not be used from sge)
test  $# -eq  0  &&   source ../D-Skel/confirm.inc

if [  $# -eq  0 ] ; then
    confirm  ../D-Skel/Smash/SkelDir
fi

sel=1

while [ $sel -ne 0 ]; do
    echo " "
    echo "select number"
    echo "0- quit"
    echo "1- edit jsfSkel"
    echo "2- submit job by bllsubmit jsfSkel"
    echo "3- see job class status"
    read sel
    if [ $sel -eq 1 ]; then
	emacs ./jsfSkel
    elif [ $sel -eq 2 ]; then
	bllsubmit jsfSkel
	echo  "See Smash/skelmemo for statsu/result"
	sleep 3
    elif [ $sel -eq 3 ]; then
	bllclass
    elif [ $sel -eq 0 ]; then
	break;
    fi
done
