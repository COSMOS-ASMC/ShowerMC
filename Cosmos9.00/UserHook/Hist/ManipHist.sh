#!/bin/sh -f
#  manipulate histograms
source $COSMOSTOP/UserHook/DisPara/Smash/setupenv.sh $0
source $COSMOSTOP/UserHook/DisPara/FleshHist/setupenv.sh $0


echo " execid =  $EXECID"
echo " fleshdir= $FLESHDIR"

###  reset next if above  is ng.
#EXECID=p3x15cos1.0
#export EXECID


source ./setupManipHistEnv.sh
source $COSMOSTOP/Scrpt/setarch.sh

for job in $JOBTYPE
do
    case $job in
	1) 
	     make -f bin2bin.mk
	    ./bin2bin$ARCH
	    ;;
	2) 
	     make -f bin2ascii.mk
	    ./bin2ascii$ARCH > $AHISTFILE
	    ;;
	    
	3|5|7)
	     make -f bin2ascii.mk
	    export -n HYBFILE0
	    echo "HYBFILE0 is not needed so it is unexported"
	    ./bin2ascii$ARCH > $AHISTFILE
	    ;;
	4|6)
	     x=`uname -o`
	     y=`echo $x | awk '{if( $1~ "Linux") print "ng"; else print "ok"}'`
	     if [ $y = "ng" ]; then
	        echo "If you are using  RedHat, FedraCore, SunOS.."
		echo "and  if histograms are too many, this process may fail"
		echo "in the mid way"
		echo "Typical error message is :fatal: can't redirect to ..."
		echo "If it happens, you may try MacOS X, Free BSD, SGI"
		echo "You need only awk script in $COSMOSTOP/Scrpt/splithisto.awk "
		echo "and ascii hist file. Usage is simple"
		echo "awk -f splithisto.awk maindir=dir ascii-histfile"
		echo "maindir is the output directory where files are goind to be stored"
	    fi
	    splithisto.sh $PLOTDIR $AHISTFILE
	    ;;
    esac
done
