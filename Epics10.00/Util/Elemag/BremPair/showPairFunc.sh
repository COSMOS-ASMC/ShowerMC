#!/bin/bash
work=/tmp/$USER/Work
if [ $#  -lt  4 ] || [ $# -gt 7 ]; then
cat <<EOF
 Usage: ./showPairFunc.sh  media  E1 E2 step [norm func LPM]
   This shows pair funciton (f=  dsigma/dk (mb) but
   the unit depend on norm )
   where k= Ek/(Eg-2me); Ek being the electron kinetic energy.
   k vs f is kept in Work/fpair.func  as numerical data
   and graph is in Work/fpair.png. 
   (old ones, if any,  will be deleted)

   media: such as Fe BGO Air ...
   E1 E2: Energy of gamma; form E1 to E2 (> 1.022e-3)
    step:  log10 step of energy
    norm:  normalization.   (default =1)
           1--> / r.l 
           2--> mb/ingredient
           3--> /(g/cm^2)
           4--> /cm
           5--> area normlization.  Integral (kmin,1) of dsigma/dk=1
    func:  Pair func by Tsai-->1; by Nelson-->2(default is 2)"
    LPM:   f no LPM, t at H.E, LPM is considered(defatul is t).
EOF
exit
fi

mkdir -p $work

lt=""   
lw=2
media=$1;  E1=$2; E2=$3; step=$4;
norm=1 
func=2 
if [ $# -ge 5 ]; then
    norm=$5
fi
if [ $# -ge 6 ]; then
    func=$6
fi
if [ $# -eq 7 ]; then
    LPM=$7
else
    LPM="t"
fi
# echo input: $media $nevent $E1 $E2 $step

nbin=`awk  'END {nbin=int(log(E2*1.0001/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 E2=$E2 step=$step /dev/null`
Eg=$E1

rm -f $work/fpair.png
rm -f $work/fpair*.func

make -f DrawPairFunc${func}.mk
if [ $? != 0 ]; then
    echo compile error
    exit 1
fi
for f in $(seq 1 $nbin); do
  echo ${f}-th energy=$Eg 
  echo $norm $media $Eg $LPM | ./drawpair${func}.out >> $work/fpair${f}.func
  echo " " >> $work/fpair${f}.func
  Eg=`awk 'END {Eg=Eg*10.**step;print Eg}' Eg=$Eg step=$step /dev/null`
done

echo "Eg=${E1}" > $work/gnuplot.com
echo "step=${step}"  >> $work/gnuplot.com
echo "media=\"${media}\""  >> $work/gnuplot.com
cat fpair0.gp >> $work/gnuplot.com

echo " plot  \"$work/fpair1.func\" w l lw ${lw} ${lt}  " >> $work/gnuplot.com
for f in $(seq 2 $nbin); do
   echo " rep  \"$work/fpair${f}.func\" w l lw ${lw} $lt " >> $work/gnuplot.com
done
cat fpair1.gp >> $work/gnuplot.com

echo  set output  \"$work/fpair.png\" >>$work/gnuplot.com
echo  rep  >>$work/gnuplot.com

gnuplot $work/gnuplot.com
