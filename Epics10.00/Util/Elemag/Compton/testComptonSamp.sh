#!/bin/bash
work=/tmp/$USER/Work
if [ $#  -ne 6 ]; then
cat <<EOF
    Usage: ./testComptonSamp.sh  media events E1 E2 step UseLogScale
      media: such as Fe BGO Air ...
     events: # of events (samplings) for each energy
      E1 E2: energy of photon; form E1 to at least E2"
       step: log10 step of energy
   UseLogScale:  0-->ordianry scale v vs ds/dv (v=Eg'/Eg)
                 1-->Log v  vs  v*ds/dv  
EOF
    exit
fi

mkdir -p work

media=$1; nevent=$2; E1=$3; E2=$4; step=$5; log=$6;

nbin=`awk  'END {nbin=int(log(E2/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 E2=$E2 step=$step /dev/null`
#echo nbin is $nbin
#exit
Ek=$E1

rm -f $work/compton*.png

make -f ComptonSamp.mk
make -f  drawCompFunc.mk
for f in $(seq 1 $nbin); do
  echo ${f}-th energy=$Ek 

  echo $Ek 1 $media  | ./drawCompFunc.out > $work/compton.func
#   get total normalization const. ( 1/(total rob/r.l))
  Nc=`awk '{print $7; exit}' $work/compton.func` 
  echo "Nc =" $Nc
  ANc=`awk '{print $8; exit}' $work/compton.func` 
  echo "ANc =" $ANc
#   output compton.func is Nc*ds/dx  ( prob/r.l)
#     dN/dx/N= ds/dx/tprob
#     trpob/N * dN/dx = ds/dx
#    Nc*tprob/N * dN/dx = Nc*ds/dx = output above 
#     ANc = tprob  

  Ncc=`echo  $ANc $Nc  $nevent | awk '{print $3/$1/$2}'`   
echo "Ncc=" $Ncc
  vmin=` awk 'END {xmin=0.9/(1.+2*Eg/0.511e-3); if(xmin > 1.e-4) xmin=1.e-4;print xmin}' Eg=$Ek /dev/null`
  echo $Ek $nevent $media   | ./a.out > $work/compton${f}.dat
  
  awk 'NR>1 {print $1}' $work/compton${f}.dat | histo -l $vmin  0.01 $Ncc  > $work/compton.hist
#  echo 5 $media $Ek $LPMeffect | ./drawbrem${func}.out > $work/brems${f}.func


  echo "Eg=${Ek}" > $work/gnuplot.com
  if [ $log -eq 0 ]; then
      echo "unset log" >> $work/gnuplot.com
      echo "set ylab \"ds/dv(/r.l)\"" >> $work/gnuplot.com
      cat compton.gp | awk '{gsub(/WWW/,this);print}' this=$work ; >> $work/gnuplot.com
  else
      echo "set log x" >> $work/gnuplot.com
      echo "set ylab \"v*ds/dv(/r.l)\"" >> $work/gnuplot.com
      echo "set format x \"%.1e\"" >> $work/gnuplot.com
      cat comptonLog.gp | awk '{gsub(/WWW/, this); print}' this=$work  >> $work/gnuplot.com
  fi
 


  echo  set output  \"$work/compton{f}.png\" >>$work/gnuplot.com
  echo  rep  >>$work/gnuplot.com

  gnuplot $work/gnuplot.com

  Ek=`awk 'END {Ek=Ek*10.**step;print Ek}' Ek=$Ek step=$step /dev/null`
done
