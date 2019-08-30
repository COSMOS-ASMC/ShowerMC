#!/bin/bash 
#  see finished job and put hostname.job# in hostnumfile
#  if cpu's < required number, exit status is 1
#  else  0
if [ $# -ne 1 ]; then
    echo "Usage: getfinishedhostnum.sh min#_of_CPU_needed"
    exit
fi
source ../Smash/setupenv.sh  $0
source ../$FLESHDIR/setupenv.sh $0
if [ -f finished ]; then
    rm -f finished
fi
for  f in $ERRDIR/*.err 
do 
  x=`tail -n 4  $f | grep "###end of run###"`
  if [ "x$x" != "x" ];  then
    echo $f >> finished
  fi
done
#
if [  -f finished ]; then
#        get hostname.0013  etc in hostnumfile
    awk -F- '{print $NF}' finished | awk -F. '{print $1"."$2}' > hostnumfile
    num=`wc -l hostnumfile | awk '{print $1}' `
    if [ $num -ge $1 ]; then
	exit 0
    else
	exit 1
    fi
else
    echo "no jobs have been completed yet"
    exit 1
fi
