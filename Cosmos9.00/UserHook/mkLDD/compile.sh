#!/bin/bash

if [ $# -ne 3 ];then
    exit
fi

srcListFile=$1
headerListFile=$2
mode=$3

WKD=`pwd`
#WKD=`pwd`/Work

CPP="ifort -DPCLinuxIFC -P -E "
#FC="ifort -g -DPCLinuxIFC -Vaxlib -align dcommons -fpconstant"
FC="ifort -DPCLinuxIFC -Vaxlib -align dcommons"
CC="icc -DPCLinuxIFC "
INC="-I${WKD}/include -I${WKD}/include/BlockData "


if [ ! -d src ];then mkdir src ; fi
if [ ! -d lib ];then mkdir lib ; fi
if [ ! -d include ];then mkdir include ; fi


if [ $mode -eq 1 ];then

    rm -f src/*.[fco]
    rm -f lib/*.[a]
    rm -f *.dyn *.dpi *.spi

    for i in `cat $headerListFile`
    do
	hfile=${COSMOSTOP}/$i
	echo -e -n "copying $hfile                            \r"
	cp -rp $hfile ${WKD}/include/
    done
    cp -fp *.h include/
    echo

    for i in `grep "\.c" $srcListFile`
      do
      cfile=${COSMOSTOP}/$i
      cp -p $cfile src/
      cd src
      $CC -prof-genx -prof-dir $WKD -c $(basename $cfile)
      cd $WKD
    done


    for i in `grep "\.f" $srcListFile`
      do
      ffile=${COSMOSTOP}/$i
      echo -e -n "preprocessing $ffile                        \r"
      cd $(dirname $ffile)
      $CPP $INC $ffile > ${WKD}/src/$(basename $ffile)
      err=$?
      cd $WKD
      if [ $err -ne 0 ];then
	  echo error
	  exit
      fi
    done
    for i in `ls ${WKD}/*.f`
    do
      ffile=${WKD}/$(basename $i)
      echo -e -n "preprocessing $ffile                          \r"
      cd $(dirname $ffile)
      $CPP $INC $ffile > ${WKD}/src/$(basename $ffile)
      err=$?
      cd $WKD
      if [ $err -ne 0 ];then
	  echo error
	  exit
      fi
    done
    echo


    echo compiling
    cd ${WKD}/src
    $FC -prof-genx -prof-dir $WKD $INC -c *.f
    mv chook.o ${WKD}
    ar r ../lib/libcosmos.a *.o
    cd $WKD

    ifort -Llib chook.o -o cosmosPCLinuxIFC -lcosmos




elif [ $mode -eq 2 ];then

    bash run_cosmos_sub.sh qgsjet2 p  1.000e+18 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh qgsjet2 p  1.000e+19 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh qgsjet2 p  1.000e+20 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh qgsjet2 Fe 1.000e+18 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh qgsjet2 Fe 1.000e+19 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh qgsjet2 Fe 1.000e+20 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh dpmjet3 p  1.000e+18 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh dpmjet3 p  1.000e+19 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh dpmjet3 p  1.000e+20 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh dpmjet3 Fe 1.000e+18 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh dpmjet3 Fe 1.000e+19 1.000 0 1.0e7 1 tmp
    bash run_cosmos_sub.sh dpmjet3 Fe 1.000e+20 1.000 0 1.0e7 1 tmp
    rm -f tmp/* fnew-* forg-* 



elif [ $mode -eq 3 ];then

    rm -f src/*.[o]

    for i in `grep "\.c" $srcListFile`
      do
      cfile=${COSMOSTOP}/$i
      cp -p $cfile src/
      cd src
      $CC -prof-use -prof-dir $WKD -c $cfile
      cd $WKD
    done

    if [ -f ${WKD}/cosmosPCLinuxIFC ];then
	mv ${WKD}/cosmosPCLinuxIFC ${WKD}/cosmosPCLinuxIFC_PGO
    fi

    cd ${WKD}/src
    $FC -prof-use -prof-dir $WKD -O3 -ipo $INC -c *.f
    $FC -prof-use -prof-dir $WKD -O3 -ipo $INC -o ${WKD}/cosmosPCLinuxIFC *.o
    cd $WKD



elif [ $mode -eq 7 ];then

    bash compile.sh $1 $2 1
    bash compile.sh $1 $2 2
    bash compile.sh $1 $2 3



elif [ $mode -eq 0 ];then

    for i in `cat $headerListFile`
    do
	hfile=${COSMOSTOP}/$i
	echo copying $hfile
	cp -rp $hfile ${WKD}/include/
    done
    cp -p *.h include/

    for i in `grep "\.c" $srcListFile`
      do
      cfile=${COSMOSTOP}/$i
      cd $(dirname $cfile)
      $CC -prof-genx -prof-dir $WKD -c $cfile
      mv `echo $cfile | sed s/"\.c"/"\.o"/g` ${WKD}/src/
      cd $WKD
    done

    for i in `grep "\.f" $srcListFile`
      do
      ffile=${COSMOSTOP}/$i
      echo preprocessing $ffile
      cd $(dirname $ffile)
      $CPP $INC $ffile > ${WKD}/src/$(basename $ffile)
      cd $WKD
    done

    cd ${WKD}/src
    $FC -O0 $INC -c *.f
    mv chook.o ${WKD}
    ar r ../lib/libcosmos.a *.o
    cd $WKD


    echo please link by yourself as follows:
    echo \$ ifort -O0 -Llib chook.o -o cosmosPCLinuxIFCO0 -lcosmos


fi

