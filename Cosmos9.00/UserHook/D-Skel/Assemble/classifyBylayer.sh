#!/bin/bash
# give arg1 /FDDbku/...
GOUTDIR=/Work1/kasahara/ForFDD
MCPU=50

exec 3<&0  <hostnumfile
nc=0



while read  host_num
do
 echo "host_num is" $host_num
#   next for is to expande *; no pluaral candidate should exist
#   in the files
#    ( from= $GOUTDIR/*$host_num".dat-r" dose not expand *)
  for to in $GOUTDIR/*$host_num".dat"  ; do
      echo "to is " $to
#  to=/XXX/yyyy/qgsjet2.tasim509.00041.dat
      basename=${to##*/}
      basename=qgsjet2p19cos0.750.dat
      echo "basename is " $basename
#  basename=qgsjet2.tasim509.00041.dat
      filen=${basename%.*}
      echo "filen is " $filen
#  filen=qgsjet2.tasim509.00041
#  mv  $GOUTDIR/*$host_num".dat-r"   $GOUTDIR/*$host_num".dat"
#      mv  $from   $to
#  $TAMCDB/bin/md2WebSysDat $GOUTDIR/  *$host_num   $1/  1
      mv $to $GOUTDIR/$basename
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



