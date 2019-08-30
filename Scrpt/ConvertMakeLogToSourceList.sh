#!/bin/bash

if [ $# -ne 1 ] || [ ! -f $1  ]
then
    exit
fi

ntot=0

dir=$(readlink -f $(dirname $0)/..)
list=`find $dir -type l`
for link in $list
do
    if [ -d $(readlink -f $link) ]
    then
	src[$ntot]=$link
	dst[$ntot]=$(readlink -f $link)
	ntot=$[$ntot+1]
    fi
done

list=`gawk '{if($2=="Entering")dir=$NF;if($1~"cppFC")print dir "/" $NF}' $1 | sed "s/'//g"`
for i in $list
do
    file=`readlink -f $i`
    out=$file
    j=0
    while [ $j -lt ${#dst[@]} ]
    do
	if [[ $file =~ "Hidden" ]] && [[ $file =~ ${dst[$j]} ]]
	then
	    out=${file/${dst[$j]}/${src[$j]}}"\tlink ${src[$j]} to\t${dst[$j]}"
	    break
	fi
	j=$[$j+1]
    done
    echo -e $out | \
	sed s@$COSMOSTOP@'${COSMOSTOP}'@g | \
	sed s@$LIBLOFT@'${LIBLOFT}'@g 
done | sort
