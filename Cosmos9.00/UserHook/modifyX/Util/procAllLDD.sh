#!/bin/bash
#
if [ $# -ne 2 ]; then
    cat <<EOF

    Usage: procAllLDD.sh   -t.hist-containing-dir dir-to-store-output

  This is to process all *-t.hist files in a given 
  directory and produce *.time files in a specified
  director.  In the input directory, *.hyb files must
  also exist.  The output directory can be the same as
  the input directory (the user must have write permission).
  
  If a file 'tList' in the output directory exists, it is
  assumed that it containes the list of files to be 
  processed; otherwise, tList is created to caintain all
  the -t.hist files.  When we finish processing of a file,
  we put "= "  at  the top of the file name; that is 
  if the first column in a tList row is "= ", the file
  in that row will not be processed

  *.time file contains fitting coefficients for (r,T10) 
  values (for all ptcl code and fai indexes).
EOF
    exit 1
fi
#  test if writable to output dir
indir=$1
outdir=$2
if [ ! -d $outdir ]; then
    echo $outdir not exist or not a directory
    exit  1
else 
    touch $outdir/dummy$$
    if [ $? -ne 0 ]; then
	echo $outdir is not writable
	exit 1
    fi
    rm -f $outdir/dummy$$
fi
#  examine append mode or not. 
#  in case of append mode, there may be newly
#  added data
pdir=`pwd`
(cd $indir; ls *-t.hist > $pdir/tListtemp )
if [ ! -s $outdir/tList  ]; then
    mv tListtemp $outdir/tList
else
    for f in `cat tListtemp`; do
	need=`grep -c $f $outdir/tList`
	if [ $need -eq  0 ]; then
	    echo $f >> $outdir/tList
	fi
    done
fi

# set arch
ARCH=`./setarch.sh`

#if [ ! -f getBasicHistoInfo${ARCH} ]; then
    make -f getBasicHistoInfo.mk
    if [ $? -ne 0 ] ; then
	echo "Unable to create  getBasicHistoInfo${ARCH}"
	exit 1
    fi
#fi
#if [ ! -f procTime${ARCH} ]; then
    make -f procTime.mk
    if [ $? -ne 0 ]; then
	echo "Unable to create procTime${ARCH} "
	exit 1
    fi
#fi


ntotal=`awk '{nl=NR}; END {print nl}' $outdir/tList`
echo "total # of -t.hist is "  $ntotal

if [ -f ./fitted.data ]; then
    rm -f ./fitted.data
fi

maxage=`awk '$1=="maxage" {print $2}' ./baseInfo`

nc=0
while [ $nc -lt $ntotal ]; do
    let nc++
    tfile=`awk 'NR==nc {print $1;exit}' nc=$nc $outdir/tList`
    if [ $tfile = "=" ]; then
	continue
    fi
    basename=`echo $tfile | awk '{i=index($1,"-t.hist"); print substr($1,1,i-1)}'`
    hybfile=${indir}/${basename}.hyb
    


    ascii=(`file -b $indir/$tfile`)
    if [ ${ascii[0]} = "ASCII" ] ; then
	format=1
	echo  "$tfile is ascii ( ${ascii[0]} )"
    else
	format=2
	echo  "$tfile is binary ( ${ascii[0]} )"
    fi
#                         0 is for -t.hist from mkLDD
    ./getBasicHistoInfo${ARCH} $format  0  $indir/$tfile  tempmemo
    if [ $? -ne 0 ]; then
	echo "fail to get basic histo info"
	exit 1
    fi

#       a bug in k90whist1.f fails to write layer number
#       correctly and 0 may be seen; we have set 30 at 
#       simulation  so  30 below is probalby ok
    layer=`awk '{if($1==0) $1=30; print $1}' tempmemo`
    depth=`awk '{print $2}' tempmemo`
    
    cosz=`awk 'NR==1 {print $9;exit}' $hybfile`
#    next two may not be needed
    erg=`awk 'NR==1 {print $6;exit}' $hybfile`      
    pcode=`awk 'NR==1 {print $3;exit}' $hybfile`      
    mu=`awk '$1=="t" && $2==layer {print $4}' layer=$layer  $hybfile`      
    echo "file#=$nc  1ry=$pcode E0=$erg cos=$cosz layer=$layer mu=$mu"

    judgeage=`awk '$1=="t" && $2==ly {if($5>=maxage) print "skip";else print "do" ;exit}' ly=$layer maxage=$maxage $hybfile`

    if [ $judgeage = "skip" ]; then
	continue
    fi

    ./timeFit.sh   ${indir}/$tfile  $cosz  $layer $mu

    if [ $? -ne 0 ]; then
	echo error in timeFit.sh
	exit 1
    fi
    mv ./fitted.data $outdir/${basename}.time
    awk 'NR != nc {print;next}; \
       NR==nc {print "= " $0}' nc=$nc  $outdir/tList > ./templist
    mv ./templist $outdir/tList
done
