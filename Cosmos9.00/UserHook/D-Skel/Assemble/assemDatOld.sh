#!/bin/bash
#  This uses $TAMCDB/bin/md2WebSysDat binary so that the  host
#  must see it and must be able to execute it.
#   Assemble ascci .dat-r in /tmp/... of various hosts
#   We classify  the data by layer number and make files
#   in the destination directory.
#   The resultant files can be put into DB by using
#   /$TAMCDB/sbin/mdCreateFDDBnew.csh
#


if [ $# -ne  1 ] ; then
    cat <<EOF
    Usage: ./assemDat.sh  destinationDir 
    destinationDir is a directory to store concatinated xxx-Lab.dat
    files.
    As xxx, EXECID is used.
    You need working directory to store the data files fetched from
    various hosts.  It must not be the same as destinationDir.
    Both directories must be able to store ~35 GB data 
EOF
    exit
fi

source ../Smash/setupenv.sh $0
#  we get OUTDIR
source ../$FLESHDIR/setupenv.sh $0
source ./setupenvHyb.sh $0
x=` echo  $OUTDIR | awk '$0 ~ "^/tmp"{print "tmp"}' `
if [ "x$x" = "xtmp" ] ; then
#     directory is /tmp/... in each host so you must gather output
#     to /tmp/... of this host
    echo "Output seems  in $OUTDIR  of each  host"
    yesno=0
    while [ x$yesno != "xy" ]
    do
      echo "You have to gather .dat-r  data into some dir in current host"
      echo "Enter  such dir (default is $OUTDIR of the current host)"
      read  dirname
      if [ x$dirname = "x" ]; then
	 GOUTDIR=$OUTDIR
	 yesno="y"
      else
         GOUTDIR=$dirname
        if [ -d $GOUTDIR ]; then 
	  echo "The dir is $GOUTDIR, is it OK ?; enter y or n"
          read yesno
        else
	  echo "$GOUTDIR not exist"
        fi
      fi
     done
  
      if [ ! -d $GOUTDIR ]; then
	 mkdir -p $GOUTDIR
	 tmp="new"
      else
	 some="`ls $GOUTDIR/`"
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
	 echo "Now we are going to gather .dat-r files in many hosts to $GOUTDIR of this host"  
	 echo "You have some files in $GOUTDIR of the current host"
	 echo "1--Delete all files in $GOUTDIR before gathering files (normal)"
	 echo "2--Delete only some files specifying file extesion(s)."
	 echo "3--Keep all files in $GOUTDIR and gather the files"
         echo "4--Files have been already gathered so keep them and proceed"
	 echo "5--Keep all files in $GOUTDIR and quit"
	 echo "Select number"
	 read num
	 test -z $num && num=0
	done
	if [ $num -eq 1 ]; then
#	    rm -f $OUTDIR/*
	    (cd $GOUTDIR; rm -f `ls`)
	elif [ $num -eq 2 ]; then
	    echo "Enter file extensions such as .hyb .hist .."
	    read extens
	    ../rmOnlySpecExt.csh $GOUTDIR "$extens"
	elif [ $num -eq 5 ]; then
	    exit
	fi
    fi
##   copy the .nrfai and .hyb files to the working directory
    if [ -f ./${EXECID}.nrfai ]; then
	cp ./${EXECID}.nrfai $GOUTDIR/
    else
	echo ".nrfai not exist in the current dir"
	exit 1
    fi
    if [ -f ./${EXECID}.hyb ]; then
	cp ./${EXECID}.hyb $GOUTDIR/
    else
	echo ".hyb not exist in the current dir"
	exit 1
    fi

    if [ $num -lt 4 ]; then
	if [ -f dummyout ]; then
	    rm -f dummyout
	fi
	if [ -f dummyerr ]; then
	    rm -f dummyerr
	fi
	echo "now gathering data; takes long time"
	echo "you can see vervose messages in dummyout and dummyerr"
	echo "By control z, and bg, send the job to back ground"

	export GOUTDIR

	./getonlyfinished.sh dat-r 1> dummyout 2> dummyerr
	 if [ $? -ne 0 ]; then
	     echo "Need to collect more data or make data by Rescue job"
	     exit 1
	 fi
    fi
else
    ./getfinishedhostnum.sh $MCPU
    if [ $? -ne  0 ]; then
	echo "# of data files is not enough"
	echo "you must wait for job completion or rescue failed jobs"
	exit 1
    fi
fi

ncpu=` ls $GOUTDIR/*.dat-r | wc -l | awk '{print $1}' `
if [ $ncpu -ne $MCPU ]; then
    echo "collected files=" $ncpu " is != " $MCPU
    echo "You may probably do once more ./assembDat.sh"
    exit 1
else
    echo $ncpu files have been collected = MCPU=$MCPU
fi


echo "Now data files in $OUTDIR are being processed; files are" 
echo "classified by the layer number"
echo "It will  take some time" 
exec 3<&0  <hostnumfile
nc=0



while read  host_num
do
#   next for is to expande *; no pluaral candidate should exist
#   in the files
#    ( from= $GOUTDIR/*$host_num".dat-r" dose not expand *)
  for from in $GOUTDIR/*$host_num".dat-r"  ; do
#  from+/XXX/yyyy/qgsjet2.tasim509.00041.dat-r
      to=${from%dat-r}dat
#  to=/XXX/yyyy/qgsjet2.tasim509.00041.dat
      basename=${to##*/}
#  basename=qgsjet2.tasim509.00041.dat
      filen=${basename%.*}
#  filen=qgsjet2.tasim509.00041
#  mv  $GOUTDIR/*$host_num".dat-r"   $GOUTDIR/*$host_num".dat"
      mv  $from   $to
#  $TAMCDB/bin/md2WebSysDat $GOUTDIR/  *$host_num   $1/  1
      $TAMCDB/bin/md2WebSysDat $GOUTDIR/  $filen   $1/  1
  done
  nc=`expr $nc + 1`
  if [ $nc = $MCPU ]; then
      break
  fi
done
# restore stdin and close 3
exec 0<&3  3<&-

if [ $MCPU -ne $nc ];  then
    echo "# of cpu's not enough; there might be running jobs"
    echo "$MCPU cpu's should exist but $nc is counted"
else
    echo "All events have been successfully assmebled to $1/"
fi





filename=${basename%.*}
echo "filename: $filename"

extension=${basename##*.}
echo "extension: $extension"

dropext=${path%gz}
echo "without ext.: $dropext"
