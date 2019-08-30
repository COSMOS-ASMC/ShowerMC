#!/bin/sh 
if [  $# != 1 ];  then
    echo 'Usage ../execflesh.sh [sge|ssh]'
    echo 'This must be invoked from FleshHist'
    exit
else
    if [ $1 != "sge" ] && [ $1 != "ssh" ]; then
	echo "$1 for the argument is invalid"
	echo "At  present ssh or sge job submission is possilbe"
	exit
    fi
fi

#echo "Are you sure that primary file at skeleton making is copied to the pwd?"
#echo "Enter n, if No"
#read yesno
#test "x$yesno" == "xn" && exit


if [ -f ./Sparam ];  then
    sparam="./Sparam"
else
    echo "There is no ./Sparam file. You have to give the parameter file"
    echo "named Sparm specified  by SkeletonFile at Skeleton making time."
    exit
fi


#echo  "Assume the file specified by SkeletonFile at Skeleton making time is"
#echo  "$sparam"
#echo  "Enter y if it is so."
#read yesno
#test "x$yesno" != "xy" && exit
##   next argument is for avoiding QA in the script
source ../setupenv.sh  
source ../Smash/setupenv.sh $0
source ./setupenv.sh  $0
#  complete name of $FLESHDIR is needed
FLESHDIR=`pwd`
export FLESHDIR
echo $FLESHDIR
source ../modparam.sh  $sparam

while [ 1 ]
do
  rescue=`pwd | awk '$1~ "Rescue" {print "Rescue";exit}'`
  if [ $NCPU = $MCPU ] ; then 
      if [ x$rescue = "xRescue"  ]; then
	  echo "This is a RESCUE job."
	  if [ $MARGIN != 0 ]; then
	      echo "MARGIN is non zero, but is neglected for rescue job"
	      MARGIN=0
	      export MARGIN
	  fi
	  if [ $ENHANCE = 1 ]; then
	      echo "If you used NCPU > MCPU at fleshing time"
	      echo "It must be the ratio NCPU/MCPU in Smach/setupenv.sh"
	  else
	      echo "ENHANCE=" $ENHANCE ". This must be NCPU/MCPU of the old run"
	      echo "Enter y if OK"
	      read yesno
	      if [ x$yesno != "xy" ]; then
		  exit
	      fi
	  fi
      else
	  if [ $ENHANCE != 1 ]; then
	      echo "ENHANCE = 1 is forced " 
	  fi
      fi
      echo '1) Do you flesh all skeletons by' "$NCPU cpus listed in $HOSTLIST"
  else
      if [ x$rescue = "xRescue"  ]; then
	  echo "for Rescue job, you have to make NCPU=MCPU"
	  exit
      fi
      ENHANCE=`awk 'END  {x=ncpu*1.0/mcpu;print x}' ncpu=$NCPU mcpu=$MCPU /dev/null`
      echo "Quasi Full M.C is assumed. ENHANCE=" $ENHANCE " will be employed"
      echo '1) Do you flesh skeletons by' " $MCPU+$MARGIN  cpus listed in $THINHOSTLIST"
  fi

  echo '2) Or specify some numbers  among them  for flesh job ?'
  echo '3) Or stop here'
  echo "Enter 1, 2 or 3"
  read  job
  if [ $job = 1 ]; then
     echo "You selected 1;  Enter y, if it is correct"
     read yesno
     test "x$yesno" != "xy"  && continue
     if [ $NCPU = $MCPU ] ; then 
	 export NOOFCPU=$NCPU
	 source  ../execflesh_all.sh  $HOSTLIST $1
     else
	 export NOOFCPU=`expr $MCPU + $MARGIN` 
	 source  ../execflesh_all.sh  $THINHOSTLIST $1
     fi
  elif [ $job = 2 ]; then
     source  ../execflesh_one.sh  $HOSTLIST $1
  elif [ $job = 3 ]; then
     exit 
  else
     continue
  fi
  break
done
