#!/bin/bash
work=Work
if [ $#  -lt  5 ] || [ $# -gt 8 ]; then
cat <<EOF
 Usage: ./IntePairFunc.sh  media outfile E1 E2 step [norm func LPM]
   This gives E vs integration of  pair funciton
   (f=  dsigma/dk (mb) but  the unit depend on norm )
   where k= Ek/(Eg-2me); Ek being the electron kinetic energy.

   media: such as Fe BGO Air ...
   outfile:  E vs total XS is put here
   E1 E2: Energy of gamma; form E1 to E2 (> 1.022e-3)
    step:  log10 step of energy
    norm:  normalization.   (default =1)
           1--> / r.l 
           2--> mb/ingredient
           3--> /(g/cm^2)
           4--> /cm
    func:  Pair func by Tsai-->1; by Nelson-->2(default is 2)
           XCOM data --> 3 
    LPM:   f--> no LPM t--> LPM at high energies (defalut is t)
EOF
exit
fi
media=$1;  
outfile=$2;
E1=$3; E2=$4; step=$5;
func=$3;

if [ $# -ge 6 ]; then
    norm=$6
else
    norm=1
fi
if [ $# -ge 7 ]; then
    func=$7
else
    func=2
fi

if [ $# -ge 8 ]; then
    LPMeffect=$8
else
    LPMeffect="t"
fi

echo "norm=" $norm "func=" $func " LPMeffect=" $LPMeffect
echo " "

echo "sleeing 3 seconds"
sleep 3


Eg=$E1
make clean
make -f IntePair${func}.mk
if [ $? != 0 ]; then
    echo compile error
    exit 1
fi

echo  $media $func  $E1 $E2 $step $norm $LPMeffect  | ./IntePair${func}.out  > $outfile
