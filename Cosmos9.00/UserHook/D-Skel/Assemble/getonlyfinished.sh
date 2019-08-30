#!/bin/bash 
if [ $# -ne 1 ]; then
    echo "Usage: getonlyfinished.sh fileextension"
    echo "fileexternsion: such as (hyb hist dat nrfai) or  all"
    exit 1
fi
source ../Smash/setupenv.sh $0
#  get  finished hostnum list in hostnumfile
./getfinishedhostnum.sh  $MCPU

count=0
for  hn  in `cat hostnumfile`
do
#   for each  hostname.numb, do; get hostname first
host=`echo $hn |  awk -F. '{print $1}'`
echo $host is being inspected
 if [ $1 =  "all" ]; then
     rsync -e ssh -avz ${host}:$OUTDIR/"*-"${hn}"*"  $GOUTDIR/
##     rsync -e ssh -avz ${host}:$OUTDIR/"*-"${hn}"*"   /Work1/temp/
 else
     echo  ${host}:$OUTDIR/"*-"${hn}.$1 
     rsync -e ssh -avz ${host}:$OUTDIR/"*-"${hn}.$1  $GOUTDIR/
##     rsync -e ssh -avz ${host}:$OUTDIR/"*-"${hn}.$1  /Work1/temp/
 fi
 count=`expr $count + 1 `
 if [ $count -eq $MCPU ]; then
     break
 fi
done


if [ $count -eq $MCPU ]; then
    echo "$1 data from $MCPU cpu's has been collected in $GOUTDIR"
    exit 0
else
    echo "$1 data from $count cpu's has been collected"
    echo "but we need more data from " `expr $MCPU - $count` " cpu's"
    exit 1
fi
