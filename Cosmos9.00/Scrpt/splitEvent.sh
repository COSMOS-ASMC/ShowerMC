#!/bin/bash
#  to split histograms data file (.ahist/.achist)  with many events.
#  
if [ $# -lt 1  ]; then
    echo "Usage $0 source-ascii-hist-file {event#1 event#2}"
    echo "Split .achist/.ahist file into files each of which"
    echo "contains only 1 event"
    echo "If event#1 is not given, all events are splitted"
    echo "The files will be stored in the same directory as"
    echo "the sourfile."
    echo "The file name will be source-file-name-1.achist etc"
    exit
fi
if [ ! -f $1 ]; then
    echo "source .ahist/.achist file:  $1 not exists"
    exit
fi

sf=$1
shift
if [ $# -eq 0 ]; then
    ev1=1
else
    ev1=$1
fi

# get last event number
last=`awk '$1=="#hist1" || $1=="#hist2" || $1=="#hist3" {ev=$2};END{print ev}' $sf`
shift
if [ $# -eq 0 ]; then
    ev2=$last
else
    ev2=$1
    test $ev2 -gt $last && ev2=$last
fi

####dir=${sf%/*}  # wrong
dir=`pwd`
#  echo ${a##*/}
#  abc.ccc.ahist
file=${sf##*/}
#  x=abc.ccc.ahist
#echo ${x##*.}
# ahist
ext=${file##*.}
#echo ${x%.*}
#abc.ccc
body=${file%.*}

echo " ev1,2: " $ev1, $ev2
echo "file: " $file
echo "sf:"$sf
echo "dir:"$dir
echo "body:"$body
echo "ext:"$ext

let n=ev1
while [ $n -le $ev2 ]; do
    awk -f $COSMOSTOP/Scrpt/splitEvent.awk evn=$n $sf > $dir/$body-$n.$ext
    n=`expr $n + 1`
done
