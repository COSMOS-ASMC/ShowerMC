#!/bin/sh 
# This is for reducing  binary .dat file at each host.
# Use this script at Assemble/ (not bg job)
# the binary program to be used for this size reduction
# is the same as the one for reduceEachSize.bin.f


#
#  This script  is to reduce the size of  binary  .dat  at each host (typcall at /tmp/$USER) 
#  and make a new .dat-r file
#  you will get a new .nrfai-r  file for each .dat.
#  you have to combine them by using assemNrfai.csh in Assemble to get
#  reduced .nrfai
#  (re-naming of .nrfai-r to .nrai is needed).
source ../FleshHist/setupenv.sh $0;  make clean; make -f reduceEachSize.binbyL.mk 
nrfaifile=`ls | grep "\.nrfai"`
# echo ".nrfai="$nrfaifile
while [ "x$nrfaifile" = "x" ]
do
    echo ".nrfai file seems missing in this dir"
    echo "Enter the name of .nrfai file "
    read nrfaifile
done

if [ -f $nrfaifile ]; then
    yesno="no"
    while [ $yesno = "no" ]
    do
      echo "Is "  $nrfaifile " the .nrfai file for the current job?"
      echo "If not, enter correct name"
      read name
      if [ "x$name" = "x" ]; then
        yesno="yes"
      else
	nrfaifile=$name
      fi
    done
else
    echo $nrfaifile  file not exist
    exit
fi




if [ -f finished ] ; then
    rm -f finished
fi

for  f  in `cat hostnumfile`
do
    echo $f
    host=`echo $f | awk -F. '{print $1}' `
#    echo "host=" $host
    numb=`echo $f | awk -F. '{print $2}' `
#    echo "numb="$numb
   ssh  $host  `pwd`/reduceEach.bin.atsitebyLSlave.sh $host $numb $nrfaifile
   echo  $f >> finished
done

