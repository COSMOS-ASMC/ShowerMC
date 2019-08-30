#!/bin/bash
if [ $# -ne 1 ]; then
    cat <<EOF
    To convert old workstation type structure construct component
    into new fortran type (eg,  a.b --> a%b )
    Usage:
    dot2perc.sh  inputFortranFile
EOF
    exit
fi
vn=`awk --version | awk '{gsub(/[\.,]/,"",$3);print $3;exit}'`
if [ $vn -lt 316 ] || [ $vn -gt 700 ]; then
  cat <<EOF
  awk verson is too old.  Please use gnu awk of version >= 3.1.6
  Some old awk's  have vn > 1000  but is older than 3.1.6
EOF
  exit
fi

file=$1
if [ "$file" == "Makfile" ]; then
    cat <<EOF
    This script should not be applied to Makefile
EOF
    exit
fi

cp $file ${file}-bkup

awk -f  $COSMOSTOP/Scrpt/dot2perc.awk $file > temp1temp1
#   sed -f  $COSMOSTOP/Scrpt/dot2perc.sed temp1temp1 > temp2temp2




same=`diff -q $file temp1temp1`
if [ -z "$same" ];  then
    echo no change
    rm temp1temp1
#    rm temp2temp2
#    rm temp3temp3
    rm ${file}-bkup
else
    mv temp1temp1 $file
fi
