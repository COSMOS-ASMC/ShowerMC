#!/bin/bash
work=/tmp/$USER/Work
if [ $#  -lt  5 ] || [ $# -gt 7  ]; then
cat <<EOF
  Usage: ./testBremSamp.sh  media events E1 E2 step [logx norm]
  media: such as Fe BGO Air  Air*0.01 ...
  events: # of events (samplings) for each energy
  E1 E2: kinetic energy of muon; form E1 to at least E2
  step:  log10 step of energy
 [logx norm]: logx:  0 (default).  x-axis is ordinary. non-zero. log10
              norm: 1 (default)  
              1--> /r.l  2 --> mb 3--> /(/(g/cm2)) 4-->/cm 5--> area
EOF
exit
fi
norms=(" " "/r.l" "(mb)" "/(g/cm2)" "/cm" " area")
       
mkdir -p $work
media=$1; nevent=$2; E1=$3; E2=$4; step=$5;

xmin=1.e-4
if [ $# -ge 6 ]; then
    logx=$6
else
    logx=0
fi
if [ $# -eq 7 ]; then
    norm=$7
else
    norm=1
fi

Ek=$E1

nbin=`awk  'END {nbin=int(log(E2*1.0001/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 E2=$E2 step=$step /dev/null`
Ek=$E1

make clean

rm -f $work/brems*.png
make -f muBrSamp.mk

make -f DrawMuBremsFunc.mk
if [ $? != 0 ]; then
    echo compile error
    exit 1
fi

for f in $(seq 1 $nbin); do
  echo ${f}-th energy=$Ek 
  echo  $norm $media  $Ek  | ./drawmubrems.out > $work/brems.func
 
#   get total normalization const. ( 1/(total rob/r.l))
  Nc=`awk '{print $4; exit}' $work/brems.func` 
  echo "Nc =" $Nc
  ANc=`awk '{print $5; exit}' $work/brems.func` 
  echo "ANc =" $ANc
#   output bremas.func is Nc*ds/dx  ( prob/r.l)
#     dN/dx/N= ds/dx/tprob
#     trpob/N * dN/dx = ds/dx
#    Nc*tprob/N * dN/dx = Nc*ds/dx = output above 
#     ANc = tprob  

  Ncc=`echo  $ANc $Nc  $nevent | awk '{print $3/$1/$2}'`   
  echo "Ncc=" $Ncc
  echo  $nevent  $Ek $media  | time  ./a.out  | awk '{print $1}' | histo  -l $xmin  0.01  $Ncc  > $work/brems.hist

  
  echo media=\"$media\" >$work/gnuplot.com
  echo Ek=${Ek}  >> $work/gnuplot.com
  echo nev=$nevent >> $work/gnuplot.com
  if [ $logx -ne 0 ]; then 
      echo "set log x" >> $work/gnuplot.com
      echo set format x '"%.2e"' >> $work/gnuplot.com
  fi
  echo set ylab \"kds/dx${norms[norm]}\"  off 1.0,7 >> $work/gnuplot.com
  cat brems.gp | awk '{gsub(/WWW/, this);print}' this=$work  >> $work/gnuplot.com

  echo  set output  \"$work/brems${f}.png\" >>$work/gnuplot.com
  echo  rep  >>$work/gnuplot.com

  gnuplot $work/gnuplot.com

  Ek=`awk 'END {Ek=Ek*10.**step;print Ek}' Ek=$Ek step=$step /dev/null`
done
