#!/bin/sh
#
if [ $# -ne 2 ]; then
    cat <<EOF

    Usage: procAll.sh  -t.hist-containing-dir dir-to-store-output

  This is to prcess all *-t.hist files in a given 
  directory and produce *.time files in a specified
  director.  In the input directory, *.hyb files must
  also exist.  The output directory can be the same as
  input directory. 
  
  If a file 'tList' in the output directory exists, it is
  assumed that it containes the list of files to be 
  processed; otherwise, tList is created to caintain all
  the -t.hist files.  If we finsh processing a file
  we put " : "  next to the file name; that is 
  if the 2nd  column in tList row is ":", the file
  in that row will not be processed.

  *.time file contans fitting coefficients for (r,t) 
  values (for all ptcl code and fai indexes).
EOF
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

if [ ! -f $outdir/tList ]; then
   ls $indir/*-t.hist > $outdir/tList
fi
 

ntotal=`awk '{nl=NR}; END {print nl}' $outdir/tList`
echo "total # of -t.hist is "  $ntotal

if [ -f $outdir/fitted.data ]; then
    rm -f $outdir/fitted.data
fi


nc=0
while [ $nc -lt $ntotal ]; do
    let nc=nc+1
    tfile=`awk 'NR==nc {print $1;exit}' nc=$nc $outdir/tList`
    fin=`awk 'NR==nc {if(NF>=2) print $1;else print " "; exit}' nc=$nc $outdir/tList`
    if [ $fin = ":" ]; then
	continue
    fi
    basename=`echo $tfile | awk '{i=index($1,"-t.hist"); print substr($1,1,i-1)'`
    hybfile=${basename}.hyb
#  get 
#      due to a bug in k90whist1.f, evid is not fully contained so
#      that layer number may         30 below is probalby ok
    
    layer=`awk '$1=="#k" {if($2==0) $2=30; print $2}' $hybfile
    cosz=`awk 'NR==1 {print $9;exit}' $hybfile`
    erg=`awk 'NR==1 {print $6;exit}' $hybfile`      
    pcode=`awk 'NR==1 {print $3;exit}' $hybfile`      
    echo $nc  $pcode $erg $cosz
    ./timeFit.sh  $pcode $erg $cosz 1>out$nc 2>err$nc
    if [ $? -ne 0 ]; then
	echo error in timeFit
	exit 1
    fi
    mv $outdir/fitted.data $outdir/${basename}.time
    awk 'NR != nc {print;next}; NR==nc {print "= " $0}' nc=$nc  $outdir/tList > $outdir/templist
    mv $outdir/templist $outdir/tList
done
