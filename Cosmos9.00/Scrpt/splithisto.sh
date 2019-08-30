#!/bin/sh -f
if [ $#  -ne  2 ] ;  then
    echo "Usage: $0  maindir ascii-histfile"
    echo " maindir is the directory where all files are to be stored"
    echo " ascii-histfile is an ascii histogram file"
    exit
fi
if [ ! -d $1 ] ; then
    mkdir -p $1
else
    some="`ls $1/`"
    if [ -n "$some" ] ; then
      echo "some files in $1"
      num=0
      while [ $num -lt  1 ] || [ $num -gt 3 ]
      do
        echo "1--remove all in $1 (normal)"
        echo "2--keep all and go"
        echo "3--keep all and quit"
        echo "Select number"
        read num
        test -z $num && num=0
      done
      case  $num  in
      1)
         (cd $1; rm -fr `ls`) ;;
#                      rm $1/*  dose no work here (why?)
      2)
        : ;;
      3)
	exit ;;
      esac
    fi
fi

awk -f $COSMOSTOP/Scrpt/splithisto.awk maindir=$1 $2
