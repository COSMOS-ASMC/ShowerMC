#!/bin/bash
work=/tmp/$USER/Work
if [ $#  -lt  5 ] || [ $# -gt 7 ]; then
    cat <<EOF
    Usage: ./testPairSamp.sh  media events E1 E2 step [LPMeffect func]
      media: such as Fe BGO Air Air*0.01  ...
      events: # of events (samplings) for ehac energy
      E1 E2: kinetic energy of electron; form E1 to at least E2
       step:  log10 step of energy
  LPMeffect: t or f for LPM effect; default is t
       func:  Pair func by Tsai-->1; by Nelson-->2(default is 2)
EOF
    exit
fi

media=$1; nevent=$2; E1=$3; E2=$4; step=$5;
if [ $# -ge 6 ]; then
    LPMeffect=$6
    if [ $# -eq 7 ]; then
	func=$7
    else
	func=2
    fi
else
    LPMeffect="t"
    func=2 
fi
# echo input: $media $nevent $E1 $E2 $step

nbin=`awk  'END {nbin=int(log(E2*1.001/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 E2=$E2 step=$step /dev/null`
#echo nbin is $nbin
#exit
Eg=$E1

rm -f $work/pair*.png
make -f PrSamp.mk
make -f DrawPairFunc${func}.mk
let dnevent=${nevent}*2
for f in $(seq 1 $nbin); do
  echo ${f}-th energy=$Eg 

  echo 1 $media $Eg $LPMeffect | ./drawpair${func}.out > $work/pair.func
#   get total normalization const. ( 1/(total rob/r.l))
  Nc=`awk '{print $4; exit}' $work/pair.func`
  echo "Nc =" $Nc
  ANc=`awk '{print $5; exit}' $work/pair.func`
  echo "ANc =" $ANc
# output pair func is Nc*ds/dx (prob/r.l)  (x=Ek/(Eg-2me))
#     dN/dx/N= ds/dx/tprob 
#     trpob/N * dN/dx = ds/dx
#    Nc*tprob/N * dN/dx = Nc*ds/dx = output above
#     ANc = tprob 
  Ncc=`echo  $ANc $Nc  $dnevent | awk '{print $3/$1/$2}'`
  echo "Ncc=" $Ncc
  
echo $nevent $Eg $media $LPMeffect  |  ./a.out  | awk '{print $2}' | histo     0. 0.01 $Ncc > $work/pair.hist


  echo "Eg=${Eg}" > $work/gnuplot.com

  cat pair.gp | awk '{gsub(/WWW/,this);print}' this=$work  >> $work/gnuplot.com

  echo  set output  \"$work/pair${f}.png\" >>$work/gnuplot.com
  echo  rep  >>$work/gnuplot.com

  gnuplot $work/gnuplot.com

  Eg=`awk 'END {Eg=Eg*10.**step;print Eg}' Eg=$Eg step=$step /dev/null`
done
