#!/bin/bash
if [ $# != 2 ]; then
    cat <<EOF
  Usage: ./getcosdist.sh input file outputfileID
EOF
    exit
fi
input=$1
output=$2
awk '($1=="izen" || $1=="ezen") && $3==1 {print $4}' $input | histo -1 0.02 > ${output}-sfn1.hist

awk '($1=="izen" || $1=="ezen") && $3==2 {print $4}' $input | histo -1 0.02 > ${output}-sfn2.hist

awk '($1=="izen" || $1=="ezen") && $3==3 {print $4}' $input | histo -1 0.02 > ${output}-sfn3.hist

awk '($1=="izen" || $1=="ezen") && $3==4 {print $4}' $input | histo -1 0.02 > ${output}-sfn4.hist

awk '($1=="izen" || $1=="ezen") && $3==5 {print $4}' $input | histo -1 0.02 > ${output}-sfn5.hist
