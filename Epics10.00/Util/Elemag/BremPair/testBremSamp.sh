#!/bin/bash
work=/tmp/$USER/Work
if [ $#  -lt  5 ] || [ $# -gt 9  ]; then
cat <<EOF
  Usage: ./testBremSamp.sh  media events E1 E2 step [LPMeffect Flpm logx func]
  media: such as Fe BGO Air  Air*0.01 ...
  events: # of events (samplings) for each energy
  E1 E2: kinetic energy of electron; form E1 to at least E2
  step:  log10 step of energy
  LPMeffect: t or f for LPM effect; default is t
  Flpm: defautl is 1.  min. is 1.  LPM is activated when Ee > 
        (minimum possilbe LPM energy) *  Flpm
  logx: 0 (default).  x-axis is ordinary. non-zero. log10
  func: one of ?, seltzer, ps1, ps3, TsaiCS to fix the
        brems function. 
        If ?, a relevant function is selected depending on
        the energy.
         
        ***Note*** 
        If you specify a function other than ?, the function 
        drawn overlayed on the  sampled histogram might be
        diffrent from the  one used in the sampling. 

        ps1--> partial screening function by Tsai
        ps3--> P.S func.  by Messel & Crowford with modification
        For each function, relevant energy region exists but
        when ? is given, it is neglected (exception: in case of 
        Seltzer function, max is  limitted to 10GeV).
EOF
exit
fi

media=$1; nevent=$2; E1=$3; E2=$4; step=$5;
mkdir -p $work
how=-1
norm=1
xmin=1.e-8
if [ $# -ge 6 ]; then
    LPMeffect=$6
else
    LPMeffect="t"
fi
if [ $# -ge 7 ]; then
    Flpm=$7
else
    Flpm=1
fi

if [ $# -ge 8 ]; then
    logx=$8
else
    logx=0
fi
if [ $# -ge 9 ]; then
    func=$9
else
    func="?"
fi


Ek=$E1
if [ $func == "ps3" ]; then
    sel=3
else
    sel=1
fi

if [ $func == "seltzer" ]; then
   E2=` echo $E2 | awk '{if($1 >= 10.) print 9.99; else print $1}' -`
fi

# echo input: $media $nevent $E1 $E2 $step


nbin=`awk  'END {nbin=int(log(E2*1.0001/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 E2=$E2 step=$step /dev/null`
#echo nbin is $nbin
#exit
Ek=$E1

make clean

rm -f $work/brems*.png
make -f BrSamp.mk

make -f DrawBremsFunc${sel}.mk
if [ $? != 0 ]; then
    echo compile error
    exit 1
fi

for f in $(seq 1 $nbin); do
  echo ${f}-th energy=$Ek 
#  echo $nevent $Ek $media $LPMeffect  |  ./a.out  | awk '{print $1}' | histo -l 1.e-5 0.01 $nevent > $work/brems${f}.hist
#  echo  $media $func $Ek $norm $xmin
  echo  $media $func $how  $Ek $norm $xmin $LPMeffect $Flpm  | ./drawbrems${sel}.out > $work/brems.func
 

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
#  echo $nevent $Ek $media $LPMeffect  |time  ./a.out  | awk '{print $1}' | histo -l 1.e-5 0.01 $Ncc  > $work/brems.hist
  echo  $nevent  $Ek $media $LPMeffect $Flpm  | time  ./a.out  | awk '{print $1}' | histo  -l $xmin  0.01  $Ncc  > $work/brems.hist



  echo "Ek=${Ek}" > $work/gnuplot.com
  if [ $logx != 0 ]; then 
      echo "set log xy" >> $work/gnuplot.com
  fi
  cat brems.gp | awk '{gsub(/WWW/, this); print}' this=$work  >> $work/gnuplot.com

  echo  set output  \"$work/brems${f}.png\" >>$work/gnuplot.com
  echo  rep  >>$work/gnuplot.com

  gnuplot $work/gnuplot.com

  Ek=`awk 'END {Ek=Ek*10.**step;print Ek}' Ek=$Ek step=$step /dev/null`
done

