#!/bin/sh
#   assemble ascci .dat-r in /tmp/... of various hosts
#   concatiante directory.
#   You may make the process as background job
#   
if [ $# -ne  1 -a  $# -ne 2 ] ; then
    echo "Usage: ./assemDatNew.sh  destinationDir [combined_filename]"
    echo "destinationDir is a directory to store concatinated .dat(=combiend_filename)"
    echo "If combined_filename is not given, EXECID is used"
    exit
fi

source ../Smash/setupenv.sh $0
#  we get OUTDIR
source ../$FLESHDIR/setupenv.sh $0
source ./setupenvHyb.sh $0
./getfinishedhostnum.sh $MCPU
if [ $? -ne  0 ]; then
  echo "# of data file is not enough"
  echo "you must wait for job completion or rescue failed jobs"
  exit 1
fi

x=` echo  $OUTDIR | awk '$0 ~ "^/tmp"{print "tmp"}' `
if [ "x$x" = "xtmp" ] ; then
#     directory is /tmp/... in each host so you must gather output
#     to /tmp/... of this host
  echo "Output seems  in $OUTDIR  of each  host"
else
  echo "Output  seems in $OUTDIR of  the current host"
  ncpu=` ls $OUTDIR/*.dat | wc -l | awk '{print $1}' `
  if [ $ncpu -ne $MCPU ]; then
      echo "collected files=" $ncpu " is != " $MCPU
     exit 1
  fi
fi
#   now hostnum$ile contains  tasim503.0350  etc

if [ $# -eq 1 ]; then
    filenm=${EXECID}.dat
else
    filenm=$2
fi
echo "Now data files in $OUTDIR are being concatinated as $1/$filenm"
echo "It will  take some time" 

if [ -f $1/$filenm ]; then
    rm -f $1/$filenm 
elif [ ! -d $1 ]; then
    echo "$1 not exist; it is created now"
    mkdir -p $1
fi



nc=0
for host_num  in `cat hostnumfile`
do
  echo "host_num= $host_num" >dummyout
 
  if [ "x$x" = "xtmp" ] ; then
     hostnm=`echo $host_num | awk -F. '{print $1}'`
     echo "data in host $hostnm is being concatinated" >> dummyout
     time ssh $hostnm cat $OUTDIR/*$host_num".dat-r" >> $1/$filenm
     retv=$?
     if [ $retv != 0 ]; then
        echo "Fatal Error: return value=" $retv >> dummyout
	exit
     fi
  else
     cat $OUTDIR/*$host_num".dat-r" >>  $1/$filenm
  fi
  nc=`expr $nc + 1`
  echo "Now data from $nc hosts has been concatinated" >> dummyout
  if [ $nc = $MCPU ]; then
      echo "nc=$nc reached MCPU=$MCPU"  >> dummyout
      break
  fi
done


if [ $MCPU -ne $nc ];  then
    echo "# of cpu's not enough; there might be running jobs" >> dummyout
    echo "$MCPU cpu's should exist but $nc is counted">> dummyout
else
    echo "All events have been successfully assmebled to $1/$filenm" >> dummyout
fi
