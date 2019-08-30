#!/bin/bash
 if [ $# != 1 ]; then
    echo "./seeFit.sh DataDir "
    exit
 fi

for f in $1/*.dat; do
    echo f is $f
    ID=`echo $f | awk '{i=index($1,"Data/"); j=index($1,"-");   \
          print substr($1,i+5,j-i-5)}'`
#    echo $ID  p or  pi or K
    
    A=`echo $f | awk '{i=index($1,"-"); j=index($1,".dat");   \
          print substr($1,i+1,j-i-1)}'`
#    echo $A     1 .. 210

    ./execFit ${ID} ${A}   < $f > temp
    cat temp >> ${ID}.coef
    cat $f > Work/xydat
    awk '{print "a"$3"="$4}' temp  > Work/gnuplot.cmd
    awk '{print "b"$3"="$5}' temp  >> Work/gnuplot.cmd
    awk '{print "c"$3"="$6}' temp  >> Work/gnuplot.cmd
    cat gnutemplate.cmd >>  Work/gnuplot.cmd

    gnuplot  Work/gnuplot.cmd

 done
