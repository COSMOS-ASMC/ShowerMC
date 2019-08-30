#!/bin/bash
maxleft=3
maxright=3
 if [ $# != 2 -a $# != 4 ]; then
    echo "./seeFit.sh hybfileDir {Ng Ne dE hNe}  [maxleft maxright] "
    exit
 elif [ $# == 4 ]; then
     maxleft=$3
     maxright=$4
 fi

for f in $1/*.hyb; do
#### file ouput;   eps file;  if pdf is usable replace to pdf
#    epsfile=${f%.hyb}.eps

    cosz=`awk '$1=="h" {print $9}' $f `
    
#    echo file is $f
    if [ $2 == "Ng" ]; then  
	awk '$1=="t" {print $3, $7}' $f  > Work/xytemp
    elif [ $2 == "Ne" ]; then  
	awk '$1=="t" {print $3, $8}' $f  > Work/xytemp
    elif [ $2 == "hNe" ]; then  
	awk '$1=="t" {print $3, $11}' $f  > Work/xytemp
    elif [ $2 == "dE" ]; then  
	awk '$1=="t" {l++; dE[l]=$13; dep[l]=$3; if(l > 1) print dep[l-1], (dE[l-1]+dE[l])/2.}' $f > Work/xytemp
    else
	echo "option " $2 " is invalid"
	exit
    fi
    ./tranFit $maxleft $maxright 2 < Work/xytemp > Work/coef

    awk '{print "a="$1}' Work/coef > Work/gnuplot.cmd
    awk '{print "b="$2}' Work/coef >> Work/gnuplot.cmd
    awk '{print "c="$3}' Work/coef >> Work/gnuplot.cmd
    awk '{print "d="$4}' Work/coef >> Work/gnuplot.cmd
    awk '{print "x0="$5}' Work/coef >> Work/gnuplot.cmd
    cat gnutemplate.cmd >>  Work/gnuplot.cmd
###  file output; change term pdf if available
#    echo "set term post eps  col " >>  Work/gnuplot.cmd
#    echo "set output \"$epsfile\"" >>  Work/gnuplot.cmd
#    echo "rep"  >>  Work/gnuplot.cmd
    gnuplot  Work/gnuplot.cmd
 done
