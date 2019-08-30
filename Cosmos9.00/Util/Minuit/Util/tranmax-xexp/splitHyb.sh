#!/bin/bash
 if [ $# != 1 ]; then
    echo "./splitHyb.sh inputfileHybFile"
    exit
 fi

 rm -f Work/file*

 awk 'BEGIN{num=0}; {print > "Work/file."num".hyb"} ; NF==0 {num++}' $1
