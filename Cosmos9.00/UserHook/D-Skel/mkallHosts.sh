#!/bin/bash
maxhosts=50000
n=1
while [ $n -le $maxhosts ]; do
    numb=`echo $n | awk 'END {printf("%5.5d\n",n)}' n=$n /dev/null`
    echo $n tasim$numb 2  1
    n=`expr $n + 1`
done
