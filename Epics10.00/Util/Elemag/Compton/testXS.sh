#!/bin/bash
work=/tmp/$USER/Work
if [ $#  -ne 3 ]; then
cat <<EOF
    Usage: ./testXS.sh  media  norm  XCOM
      media: such as Fe BGO Air ...
       norm: 1--> /r.l   2-->/(g/cm2)  3-->/cm  4 mb.
       XCOM: 0-->use ordinary Xsection
             1-->use XCOM xsection 
EOF
    exit
fi

mkdir -p $work

media=$1;  norm=$2
# from 10keV to 100GeV log step 0.1
#E1=10.e-6
#E2=100.
#step=0.1
#nbin=`awk  'END {nbin=int(log(E2/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 E2=$E2 step=$step /dev/null`
#echo nbin is $nbin
#exit
#  from 10 keV
#Ek=$E1
XCOM=$3
rm -f $work/compton*.png

make -f DrawComptonXS.mk
#for f in $(seq 1 $nbin); do
#  echo ${f}-th energy=$Ek 

  echo  $XCOM $media  | ./XS.out > $work/compton.xs

  echo "set xlab \"Eg(MeV)\"" > $work/gnuplot.com
  echo "set log xy"  >> $work/gnuplot.com
  if [ $norm -eq 1 ]; then
      ylab="/r.l"
  elif [ $norm -eq 2 ]; then
      ylab="/(g/cm2)"
  elif [ $norm -eq 3 ]; then
      ylab="/cm"
  elif [ $norm -eq 4 ]; then
      ylab="mb"
  else
      echo "error norm=",$norm
      exit 1
  fi
  echo "set ylab \"$ylab\"" >> $work/gnuplot.com

  let u=$norm+1
  echo "plot  \"$work/compton.xs\" u 1:$u  w l lw 3 lt 3" >> $work/gnuplot.com
  cat XS.gp | awk '{gsub(/WWW/,this);print}' this=$work  >> $work/gnuplot.com


  echo  set output  \"$work/compton{f}.png\" >>$work/gnuplot.com
  echo  rep  >>$work/gnuplot.com

  gnuplot $work/gnuplot.com

  Ek=`awk 'END {Ek=Ek*10.**step;print Ek}' Ek=$Ek step=$step /dev/null`

