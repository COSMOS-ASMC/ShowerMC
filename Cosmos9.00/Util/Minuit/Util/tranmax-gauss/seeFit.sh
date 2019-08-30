#!/bin/bash
 if [ $# != 1 ]; then
    echo "./seeFit.sh inputfileDir "
    exit
 fi

for f in $1/*.hyb; do
    cosz=`awk '$1=="h" {print $9}' $f `
#    echo file is $f
    awk '$1=="t" {print $3, $8}' $f  > .xytemp
#    echo "(x,y)="
    ./tranFit 2 < .xytemp > Work/coef
#    echo "fitted coef"
    awk '{print "a="$1}' Work/coef > Work/gnuplot.cmd
    awk '{print "s="$2}' Work/coef >> Work/gnuplot.cmd
    awk '{print "x0="$3}' Work/coef >> Work/gnuplot.cmd
    echo "set xr[x0-100:x0+100]"  >>Work/gnuplot.cmd
    dqm='"'
    echo "plot $dqm.xytemp$dqm ps 2 pt 7" >> Work/gnuplot.cmd
    cat gnutemplate.cmd >>  Work/gnuplot.cmd
    gnuplot  Work/gnuplot.cmd
 done
