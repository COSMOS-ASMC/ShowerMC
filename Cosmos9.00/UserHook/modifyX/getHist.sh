#!/bin/bash
if [ $# != 2 ]; then
    cat <<EOF
     Usage: ./getHist.sh 'files' id
     files such as p17CM-*.dat.  must be enclosed by ' '
     id    such as p17CM
EOF
    exit
fi

awk '$1=="xd" && $2==4 {print $4}' $1  | histo 0 0.01 20000 > ${2}-piall-x-ord.hist

awk '$1=="xd" && $2==4 {print $4}' $1  | histo -l 1.e-4 0.05 20000 > ${2}-piall-x.hist

awk '$1=="xd" && $2==4 {print $5}' $1  | histo -1 0.25 20000 > ${2}-piall-y.hist

awk '$1=="xd" && $2==4 {print $6}' $1  | histo -1 0.25 20000 > ${2}-piall-eta.hist

