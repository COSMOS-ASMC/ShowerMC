#!/bin/bash
if [ $# -ne 6 ]; then
    cat<<EOF
    Usage: getHistPxPy.sh code chg1 chg2  pz1 pz2 file
   code:  ptcl code 
   chg1,chg2:   chg1<=charge <=chg2 is selected
   pz1,pz2:     pz1< pz < pz2 is selected
   and make histogram of px and py
EOF
exit
fi
code=$1
chg1=$2
chg2=$3
pz1=$4
pz2=$5
file=$6
awk '$1==code && $3>=chg1 && $3<=chg2 && $7>pz1  && $7<pz2 \
 {print $5}'  code=$code chg1=$chg1 chg2=$chg2 pz1=$pz1 pz2=$pz2 \
 $file | histo -1.5 0.02  > $file-$code-chg$chg1-$chg2-Pz$pz1-$pz2-px.hist

awk '$1==code && $3>=chg1 && $3<=chg2 && $7>pz1  && $7<pz2 \
 {print $6}'  code=$code chg1=$chg1 chg2=$chg2 pz1=$pz1 pz2=$pz2 \
 $file | histo -1.5 0.02  > $file-$code-chg$chg1-$chg2-Pz$pz1-$pz2-py.hist


