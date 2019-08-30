#!/bin/sh
#   assemble ascci histogram  data to get unified histogram data
#

# we get NCP and HOSTLIST
source ../Smash/setupenv.sh $0
#  we get OUTDIR

#source ../$FLESHER/setupenv.sh $0
source ../$FLESHDIR/setupenv.sh $0
source ./setupenvHist.sh $0
test  -f $HISTFILET  &&    rm -f $HISTFILET 
make clean
make -f add2hist.mk
source $COSMOSTOP/Scrpt/setarch.sh
if [ ! -f  add2hist$ARCH ]; then
 echo  add2hist$ARCH not exists
 exit 1
fi
export GOUTDIR=$OUTDIR
x=` echo  $OUTDIR | awk '$0 ~ "^/tmp"{print "tmp"}' `
if [ "x$x" = "xtmp" ] ; then
#     directory is /tmp/... in each host so you must gather output
#     to /tmp/... of this host
    echo "Output seems  in $OUTDIR  of each  host"
    echo "You have to gather .hist data into $OUTDIR of this host"
    if [ ! -d $OUTDIR ]; then
	mkdir -p $OUTDIR
	tmp="new"
    else
	some="`ls $OUTDIR/`"
	if [ -n "$some" ] ; then
	    tmp="old"
	else
	    tmp="new"
	fi
    fi
    if [ $tmp = "new"  ]; then
	num=1
    else
	num=0
	while [ $num -lt 1 ] || [ $num -gt 5 ]
	do
	 echo "Now we are going to gather .hist files in many hosts to $OUTDIR of this host"  
	 echo "You have some files in $OUTDIR of the current host"
	 echo "1--Delete all files in $OUTDIR before gathering files (normal)"
	 echo "2--Delete only some files specifying file extesion(s)."
	 echo "3--Keep all files in $OUTDIR and gather the files"
         echo "4--Files have been already gathered so keep them and proceed"
	 echo "5--Keep all files in $OUTDIR and quit"
	 echo "Select number"
	 read num
	 test -z $num && num=0
	done
	if [ $num -eq 1 ]; then
#	    rm -f $OUTDIR/*
	    (cd $OUTDIR; rm -f `ls`)
	elif [ $num -eq 2 ]; then
	    echo "Enter file extensions such as .hist .hyb .."
	    read extens
	    ../rmOnlySpecExt.csh $OUTDIR "$extens"
	elif [ $num -eq 5 ]; then
	    exit
	fi
    fi
    if [ $num -lt 4 ]; then
	./getonlyfinished.sh hist
	if [ $? -ne 0 ]; then
	    echo "Need to collect more data or make data by Rescue job"
	    exit 1
	fi
    fi
else
    ./getfinishedhostnum.sh  $MCPU
    if [ $? -ne  0 ]; then
	echo "# of data file is not enough"
	echo "you must wait for job completion or rescue failed jobs"
	exit 1
    fi
fi
# chang stdin
#exec 3<&0  <$HOSTLIST
exec 3<&0  <hostnumfile


nc=0
while read host_num
do
  nc=`expr $nc + 1`
  if [ $nc -eq 1 ];  then
      cp $OUTDIR/*$host_num".hist"  $HISTFILE0
      echo "first file is copied to $HISTFILE0"
      if [ $nc -eq $MCPU ]; then
	  break
      fi
      continue
  fi
  numb=` echo $host_num | awk  -F. '{print $2}'`
  file=$numb".hist"
  echo "proceesing file with number: $numb"
  pwd=`pwd`
  cd $OUTDIR
  for f in `ls ./`
  do  
     x=${f%%$file}
     test $x = $f && continue
     break
  done
  if [ $x = $f ]; then
      echo "strange;  file with " $numb
      exit 1
  fi
  cd $pwd
  HISTFILEX=$OUTDIR/$x$file      
  export  HISTFILEX
  ./add2hist$ARCH
  mv -f $HISTFILET $HISTFILE0
  if [ $nc -eq $MCPU ]; then
      break
  fi
done

# restore stdin and close 3
exec 0<&3  3<&-

if [ $MCPU -ne $nc ];  then
   echo "# of cpu's not enough; there might be running jobs/or dead cpu"
   echo "$MCPU cpu's should exist but $nc is counted"
else
    echo "All events have been successfully assmebled to $HISTFILE0"
fi
