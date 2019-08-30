#!/bin/sh 
## don't touch next line; test is needed for sge job  (may not be used from sge)
test  $# -eq  0  &&     source ../confirm.inc

#  don't touch below.
#  if used from script, we skip the next lines.
if [  $# -eq  0 ] ; then
    confirm  ../FleshHist/SmSkelDir
fi

sel=1
while [ $sel -ne 0 ]; do
    echo " --- smash --- "
    echo "select number"
    echo "0- quit"
    echo "   You have to edit setupenvcls.sh to fix MCPU, NCPU even mpi job"
    echo "   Note MaxCpu must be  >= NCPU in Zprivate2.h"
    echo "1- look at skelmemo (by less)"
    echo "   Find a good event with moderate first collision depth"
    echo "   (in less mode, find rank by /rank  and memorize it's number)"
    echo "2- give rank number and select that event"
    echo "3- edit jsfSmash (which is a submit script)"
    echo "4- edit setupenvcls.sh"  
    echo "5- setupenvcls.sh is OK so goto ../ and adjust Hosts/thinhosts"
    echo "6- make smash executable(make -f smashSkel.mk)"
    echo "7- submit job by bllsubmit jsfSmash"
    echo "8- see job class status (by bllclass)"
    read sel
    if [ $sel -eq 1 ]; then
	less ./skelmemo
    elif [ $sel -eq 2 ]; then
	echo "Enter rank #"
	read num
	SkelFile=`awk 'END{printf("Skeleton-%d",num)}' num=$num /dev/null`
	if [ -f ./SkelDir/$SkelFile ]; then
	    cp ./SkelDir/$SkelFile ./Skeleton
	    echo "./SkelDir/$SkelFile has been copied to ./Skeleton"
	    sleep 3
	else
	    echo "There is no ./SkelDir/$SkelFile "
	    sleep 5
	fi
    elif [ $sel -eq 3 ]; then
	emacs ./jsfSmash
    elif [ $sel -eq 4 ]; then
	emacs ./setupenvcls.sh
    elif [ $sel -eq 5 ]; then
      echo "if # of ranks are large, this may take some time"
      (cd ../; ./mkHosts.csh > Hosts; ./mkThinHosts.sh)
      echo "hosts business ended"
      sleep 3
    elif [ $sel -eq 6 ]; then
	make clean;make -f smashSkel.mk
    elif [ $sel -eq 7 ]; then
	bllsubmit jsfSmash
	echo "./err shows job status"
	sleep 3
    elif [ $sel -eq 8 ]; then
	bllclass
    elif [ $sel -eq 0 ]; then
	break;
    fi
done
