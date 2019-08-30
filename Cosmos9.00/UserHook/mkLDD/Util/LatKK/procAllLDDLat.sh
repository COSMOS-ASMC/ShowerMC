#!/bin/bash
#
if [ $# -ne 3 ]; then
    cat <<EOF
    Usage: procAllLDDLat.sh  -r.hist-containing-dir dir-to-store-output  putage

  This is to process all *-r.hist files in a given 
  directory and produce *.lat files in a specified
  director.  In the input directory, *.hyb files must
  also exist.  The output directory can be the same as
  the input directory (the user must have write permission).
  putage: 0--> no age and cogd are put in the final output.
          1--> age and cogd are put in the final output.
  If a file 'rList' in the output directory exists, it is
  assumed that it containes the list of files to be 
  processed; otherwise, rList is created to caintain all
  the -r.hist file names.  When we finish processing of a file,
  we put "= "  at  the top of the file name; that is 
  if the first column in an rList row is "= ", the file
  in that row will not be processed

  *.lat file contains fitting coefficients for
   r<0.1<r<1<r<10<r<100.
    values (for all ptcl code and fai indexes).
EOF
    exit 1
fi
#  test if writable to output dir
indir=$1
outdir=$2
putage=$3
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

if [ ! -s $outdir/rList  ]; then
    pdir=`pwd`
    (cd $indir; ls *-r.hist > $pdir/rList )
    mv rList $outdir/
fi

# set arch
ARCH=`../setarch.sh`

#if [ ! -f getBasicHistoInfo${ARCH} ]; then
 make -f getBasicHistoInfo.mk
    if [ $? -ne 0 ] ; then
	echo "Unable to create  getBasicHistoInfo${ARCH}"
	exit 1
    fi
#fi
#if [ ! -f procTime${ARCH} ]; then
    make -f procLat.mk
    if [ $? -ne 0 ]; then
	echo "Unable to create procLat${ARCH} "
	exit 1
    fi
#fi


ntotal=`awk '{nl=NR}; END {print nl}' $outdir/rList`
echo "total # of -r.hist is "  $ntotal

if [ -f ./fitted.data ]; then
    rm -f ./fitted.data
fi

# maxage=`awk '$1=="maxage" {print $2}' ./baseInfo`

nc=0
while [ $nc -lt $ntotal ]; do
    let nc++
    rfile=`awk 'NR==nc {print $1;exit}' nc=$nc $outdir/rList`
    if [ $rfile = "=" ]; then
	continue
    fi
    basename=`echo $rfile | awk '{i=index($1,"-r.hist"); print substr($1,1,i-1)}'`
    hybfile=${indir}/${basename}.hyb
    

    ascii=(`file -b $indir/$rfile`)
    if [ ${ascii[0]} = "ASCII" ] ; then
	format=1
	echo  "$rfile is ascii ( ${ascii[0]} )"
    else
	format=2
	echo  "$rfile is binary ( ${ascii[0]} )"
    fi
#                         
    ./getBasicHistoInfo${ARCH} $format   $indir/$rfile  tempmemo
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
    if [ $putage -eq  1 ]; then
	age=`awk '$1=="t" && $2==ly {print $5;exit}' ly=$layer  $hybfile`
	cogd=`awk '$1=="t" && $2==ly {print $6;exit}' ly=$layer  $hybfile`
	echo age=$age, cogd=$cogd
    fi

#############
#    judgeage=`awk '$1=="t" && $2==ly {if($5>=maxage) print "skip";else print "do" ;exit}' ly=$layer maxage=$maxage $hybfile`

#    if [ $judgeage = "skip" ]; then
#	continue
#    fi
#############
    if [ $putage -eq 0 ]; then
	./latFit.sh   ${indir}/$rfile   $layer 
    else
	./latFit.sh   ${indir}/$rfile   $layer $age $cogd
    fi

    if [ $? -ne 0 ]; then
	echo error in latFit.sh
	exit 1
    fi
    mv ./fitted.data $outdir/${basename}.lat
    awk 'NR != nc {print;next}; \
       NR==nc {print "= " $0}' nc=$nc  $outdir/rList > ./templist$$
    mv ./templist$$ $outdir/rList
done
