#!/bin/bash
work=/tmp/$USER/Work
if [ $#  -lt  6 ] || [ $# -gt 13 ]; then
cat <<EOF
  This draws a specified brems function in a given
  energy region even if it is outsied of the proper
  region. 

 Usage: ./showBremFunc.sh  media func how E1 E2 step [norm
           LPMeffect Flpm logx xmin linew linet ]

   This shows brems funciton (f= k dsigma/dk (mb) but
   the unit depend on norm )
   where k= Eg/Ek; Ek being the electron kinetic energy.
   k vs f is kept in Work/fbrems.func  as numerical data
   and graph is in Work/fbrems.png. 
   (old ones, if any,  will be deleted)

   media: such as Fe BGO Air ...
   func: one of ?, seltzer, ps1, ps3, TsaiCS  to fix the
          brems function. 
          If ?, a relevant function is selected depending on
          the energy.
          ps1--> partial screening function by Tsai
          ps3--> P.S func.  by Messel & Crowford with modification
                      by Koh & Motz
          When a specific function is selected,  calculation is tried
          even if the energy is outside of the relevant energy region.
          (Exception: in case of  Seltzer function, max is  limitted 
          to 10GeV).
   how:   1, normalize from lower E cross-section
          0, no normalization
         -1, normalzie from higher E cross-section 
   E1,E2:  kinetic energy of electron, form E1 to E2 
    step:  with log10 step of energy. 
    norm: normalization.   (default =1)
           1--> / r.l 
           2--> mb/ingredient
           3--> /(g/cm^2)
           4--> /cm
           5--> area normlization.  Integral (kmin,1) of dsigma/dk=1
   LPMeffect: (default is t). If f, LPM is neglected
   Flpm:    (default is 1). (must be >=1).  LPM works when Ee
        > Flpm*(minimum LPM energy) and LPMeffect=t
   logx:   default 0.  (no log x for display). 1=> log x is used.
    xmin:  if given (>0), default value is overtaken. If 0 or not
           given, a default value for each function will be used
           Note: normalization depends on xmin.
   linew:  plot line width. default is 2.
   linet:  line type. defautl is from changing values ( 1 2...)
           -1-->all black. 0-->all black dot 1-->all red  ....  
EOF
exit
fi
norms=(" " "/r.l"  "(mb)" "/(g/cm2)" "/cm" " area")

mkdir -p $work
media=$1;  
func=$2;
how=$3
echo "media=" $media " func=" $func "how=" $how
E1=$4; E2=$5; step=$6;
echo "E1=" $E1 " E2=" $E2 " step=" $step


if [ $# -ge 7 ]; then
    norm=$7
else
    norm=1 
fi


if [ $# -ge 8 ]; then
    LPMeffect=$8
else
    LPMeffect="t"
fi

if [ $# -ge 9 ]; then
    Flpm=$9
else
    Flpm=1
fi
shift
if [ $# -ge 9 ]; then
    logx=$9
else
    logx=0
fi

shift
if [ $# -ge 9 ]; then
    xmin=$9
else
    xmin=0.  # non zero default will be used
fi

shift
if [ $# -ge 9 ]; then
    lw=$9
else
    lw=2
fi

shift
lt="" 
if [ $# -ge 9 ]; then
   lt=" lt $9"
   echo "lt:" $lt
fi  

echo "norm=" $norm "LPMeffect=" $LPMeffect "Flpm=" $Flpm  " logx=" $logx " xmin=" $xmin "lw=" $lw " lt=" ${lt} 

echo " "

echo "sleeing 3 seconds"
sleep 3

 nbin=`awk 'END {nbin=int(log(E2*1.0001/E1)/log(10.)/step)+1;print nbin}'  E1=$E1 E2=$E2 step=$step /dev/null`

if [ $func == "seltzer" ]; then
    nbin2=`awk  'END {nbin=int(log(9.999/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 step=$step /dev/null `
   if [ $nbin2 -lt $nbin ]; then
      nbin=$nbin2
   fi
fi

echo "nbin=" $nbin

Ek=$E1
if [ $func == "ps3" ]; then
    sel=3
else
    sel=1
fi
rm -f $work/fbrems*.png
rm -f $work/fbrems*.func

make clean

make -f DrawBremsFunc${sel}.mk
if [ $? != 0 ]; then
    echo compile error
    exit 1
fi

for f in $(seq 1 $nbin); do
  echo ${f}-th energy=$Ek 
  echo  "media  func how  Ek norm xmin LPMeffect Flpm:"
  echo  $media $func $how  $Ek $norm $xmin $LPMeffect $Flpm
  echo  $media $func $how  $Ek $norm $xmin $LPMeffect $Flpm | ./drawbrems${sel}.out > $work/fbrems${f}.func
  echo " " >> $work/fbrems${f}.func
  Ek=`awk 'END {Ek=Ek*10.**step;print Ek}' Ek=$Ek step=$step /dev/null`
done

echo "Ek=${E1}" > $work/gnuplot.com
echo "step=${step}"  >> $work/gnuplot.com
echo "media=\"${media}\""  >> $work/gnuplot.com

if [ ${logx} !=  0 ]; then
    echo "set log x"  >> $work/gnuplot.com
fi
 echo set ylab \"kds/dx${norms[norm]}\"  off 1.0,7 >> $work/gnuplot.com
cat fbrems0.gp  >> $work/gnuplot.com
echo " plot  \"$work/fbrems1.func\" w l lw $lw $lt " >> $work/gnuplot.com
for f in $(seq 2 $nbin); do
   echo " rep  \"$work/fbrems${f}.func\" w l lw $lw $lt " >> $work/gnuplot.com
done
cat fbrems1.gp >> $work/gnuplot.com
echo  set output  \"$work/fbrems.png\" >>$work/gnuplot.com
echo  rep  >>$work/gnuplot.com

gnuplot $work/gnuplot.com
