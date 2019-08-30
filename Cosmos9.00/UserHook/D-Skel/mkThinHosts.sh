#!/bin/bash
#  0) allHosts must exist
#  1)  from Hosts; make a table which contains random number in the last 
#   item
#  2)  sort it;  key= random number
#  3)  select needed number of hosts
if [ -f temphost ]; then
    rm -f temphost
fi
awk 'BEGIN{ srand()}; {$4=rand();print }' Hosts > temphost
if [ -f ThinHosts ]; then
    rm -f ThinHosts
fi

sort -k 4 -g temphost > ThinHosts
rm -f temphost
mv ThinHosts temphost

touch ThinHosts
source Smash/setupenvcls.sh  $0
number=$[ $MCPU + $MARGIN ]
head -n $number temphost > ThinHosts
ncount=`wc -l ThinHosts | awk '{print $1}' `
if [ $ncount -lt $number ]; then
    echo "number of cpu's is =" $ncount " < " $number
    echo "i.e, number of real cpu's is not enough"
    exit
else
    # sort by the number in the first field
    sort -n -k 1 ThinHosts > temphost
    mv temphost ThinHosts
    echo "ThinHosts has been created with " $ncount " cpu's"
fi
#  rm -f temphost
