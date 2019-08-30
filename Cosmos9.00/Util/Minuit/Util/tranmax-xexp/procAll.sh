#!/bin/bash
maxleft=3
maxright=3
 if [ $# != 2 -a  $# != 4  ]; then
cat <<EOF
Usage:   ./procAll.sh inputfileDir outputFile [maxleft maxright]
    maxleft: howmany point do you use before apparent max
    maxrigth: //                     after //
    in total maxleft+maxright+1  points will be used.
    default is 3+3+1
EOF
exit
 elif [ $# == 4 ]; then
     maxleft=$3
     maxright=$4
 fi
 rm -f $2
 for f in $1/*.hyb; do
    cosz=`awk '$1=="h" {print $9;exit}' $f `
#          gamma
    awk '$1=="t" {print $3, $7}' $f | ./tranFit $maxleft $maxright 5  > temp
    gxmax=`awk  '{print $1;exit}' temp `
    gaxmax=`awk  '{print $2;exit}' temp `
#          electron
    awk '$1=="t" {print $3, $8}' $f | ./tranFit $maxleft $maxright 5  > temp
    exmax=`awk  '{print $1;exit}' temp `
    eaxmax=`awk  '{print $2;exit}' temp `
#        hybrid  electron
    awk '$1=="t" {print $3, $11}' $f | ./tranFit $maxleft $maxright 5  > temp
    hxmax=`awk  '{print $1;exit}' temp `
    haxmax=`awk  '{print $2;exit}' temp `

#     dE  for layer l, get average at (l,l+1); cause dE at l is dE in l-1 to l.
#     for the last layer no output
#          1  
#          2      depth 1 (1+2)/2
#          3      depth 2 (2+3)/2
#          4      depth 3 (3+4)/2
#      .....
#          L-1   
#          L      depth L-1 (L-1+L)/2
    awk '$1=="t" {l++; dE[l]=$13; dep[l]=$3; if(l > 1) print dep[l-1], (dE[l-1]+dE[l])/2.}' $f \
      |	 ./tranFit $maxleft $maxright 5  > temp
    Exmax=`awk  '{print $1;exit}' temp `
    Eaxmax=`awk  '{print $2;exit}' temp `

    echo $gxmax $gaxmax $exmax $eaxmax $Exmax $Eaxmax  $hxmax $haxmax \
	|  awk '{print $1/cosz, $2/cosz, $3/cosz, $4/cosz, $5/cosz, $6/cosz, $7/cosz, $8/cosz,  cosz}' cosz=$cosz  >> $2
 done
