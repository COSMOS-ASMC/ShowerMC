#!/bin/csh
#  This makes r vs T for a given fai, code and percentage
#  Difference from mkRvsT.csh is that treats only one combination
# 
if( $#argv != 5 ) then
    echo "Usage: mkRvsT.csh fai code inputTFdata percent ouputfile"
    echo " where fai:1->(-15,15) deg. 2->(15,45)... 12->(315,345)"
    echo " code, 1,2,3,4 ->g,e,mu,h"
    echo " inputTFdata is such as tf.data"
    echo " percent 2->5%, 3->10%, 4->20%, 5->30%, 6->40%, 7->50%.."
    echo " outputfile which will contain r(m.u) and time(ns)"
    exit
endif
set head=$COSMOSTOP/UserHook/mkLDD/Util
 awk -f $head/getTimeAtFai.awk fai=$1 code=$2 $3 |  awk -f $head/getTimeAtFrac.awk nth=$4 > $5
