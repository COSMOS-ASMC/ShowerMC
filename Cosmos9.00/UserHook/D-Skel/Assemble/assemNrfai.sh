#!/bin/sh
#   assemble ascci nrfai data to get unified nrfai data
#
num=0
while [ $num -lt 1 ] || [ $num -gt 2 ]
  do
  echo "1)  now you are combining .nrfai for the first time; each .nrfai was created at each host"
  echo "or"
  echo "2)  now you are combining .nrfai-r which  has been created when reducing .dat from each host"
  echo "Select number above"
  read num
  test -z $num && num=0
done
echo "your input is" $num
if [ $num -eq 1 ]; then
    GETALL="yes"
    ext="nrfai"
else 
    GETALL="no"
    ext="nrfai-r"
fi
export GETALL

# we get NCPU and HOSTLIST etc
source ../Smash/setupenv.sh $0
if [ $NCPU = $MCPU ]; then
    XCPU=$NCPU
    XHOST=$HOSTLIST
else
    XCPU=$MCPU
    XHOST=$THINHOSTLIST
fi
export XCPU
export XHOST

#  we get OUTDIR

#source ../$FLESHER/setupenv.sh $0
source ../$FLESHDIR/setupenv.sh $0
source ./setupenvNrfai.sh $0
test  -f $NRFAIFILET  &&    rm -f $NRFAIFILET 
make clean
make -f add2nrfai.mk
source $COSMOSTOP/Scrpt/setarch.sh
if [ ! -f  add2nrfai$ARCH ]; then
 echo  add2nrfai$ARCH not exists
 exit 1
fi

export GOUTDIR=$OUTDIR

x=` echo  $OUTDIR | awk '$0 ~ "^/tmp"{print "tmp"}' `
if [ "x$x" = "xtmp" ] ; then
#     directory is /tmp/... in each host so you must gather output
#     to /tmp/... of this host
    echo "Output seems  in $OUTDIR  of each  host"
    echo "You have to gather .$ext data into $OUTDIR of this host"
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
	 echo "Now we are going to gather .$ext files in many hosts to $OUTDIR of this host"  
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
	    echo "Enter file extensions such as .nrfai .hist .."
	    read extens
	    ../rmOnlySpecExt.csh $OUTDIR "$extens"
	elif [ $num -eq 5 ]; then
	    exit
	fi
    fi
    if [ $num -lt 4 ]; then
	./getonlyfinished.sh $ext
	if [ $? -ne 0 ]; then
	    echo "Need to collect more data or make data by Rescue job"
	    exit 1
	fi
	if [ $ext = "nrfai-r" ]; then 
	    for f in $OUTDIR/*.$ext
	    do
	      mv $f ${f%nrfai-r}nrfai 
	    done
	fi
    fi
else
    ./getfinishedhostnum.sh $MCPU
    if [ $? -ne  0 ]; then
	echo "# of data file is not enough"
	echo "you must wait for job completeion or Rescue failed jobs"
	exit 1
    fi
fi

# chang stdin

#exec 3<&0  <$HOSTLIST
exec 3<&0  < hostnumfile

nc=0
while read host_num
do
  nc=`expr $nc + 1`
  if [ $nc -eq 1 ];  then
      cp $OUTDIR/*$host_num".nrfai"  $NRFAIFILE0
      echo "first file is copied to $NRFAIFILE0"
      if [ $nc -eq $MCPU ]; then
	  break
      fi
      continue
  fi
  numb=` echo $host_num | awk -F.  '{print $2}' `
  file=$numb".nrfai"
  echo "proceesing file with number: $numb"
  pwd=` pwd `
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
  NRFAIFILEX=$OUTDIR/$x$file      
  export  NRFAIFILEX
  ./add2nrfai$ARCH
  mv -f $NRFAIFILET $NRFAIFILE0
  if [ $nc -eq $MCPU ]; then
      break
  fi
done

# restore stdin and close 3
exec 0<&3  3<&-

if [ $MCPU -ne  $nc ];  then
   echo "# of cpu's not enough; there might be running jobs"
   echo "$MCPU cpu's should exist but $nc is counted"
else
    echo "All events have been successfully assmebled to $NRFAIFILE0"
fi
