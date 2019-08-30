#!/bin/bash
work=Work
if [ $#  -lt  5 ] || [ $# -gt 7 ]; then
cat <<EOF
  This draws a specified brems function in a given
  energy region even if it is outsied of the proper
  region. 

 Usage: ./showBremFunc.sh  media func E1 E2 step [norm xmin]

   This shows brems funciton (f= k dsigma/dk (mb) but
   the unit depend on norm )
   where k= Eg/Ek; Ek being the electron kinetic energy.
   k vs f is kept in Work/fbrems.func  as numerical data
   and graph is in Work/fbrems.png. 
   (old ones, if any,  will be deleted)

   media: such as Fe BGO Air ...
   func: one of ?, seltzer, ps1, ps3, TsaiCS, CS+LPM to fix the
          brems function. 
          If ?, a relevant function is selected depending on
          the energy.
          ps1--> partial screening function by Tsai
          ps3--> P.S func.  by Messel & Crowford with modification
          For each function, relevant energy region exists but
          when ? is given, it is neglected (exception: in case of 
          Seltzer function, max is  limitted to 10GeV).
   E1 E2: kinetic energy of electron; form E1 to E2
    step: log10 step of energy
    norm: normalization.   (default =1)
           1--> / r.l 
           2--> mb/ingredient
           3--> /(g/cm^2)
           4--> /cm
           5--> area normlization.  Integral (kmin,1) of dsigma/dk=1
    xmin:  if given, default value is overtaken.
           Note: normalization depends on xmin.
EOF
exit
fi
   
media=$1;  
func=$2;
shift; 
E1=$2; E2=$3; step=$4;
norm=1 
xmin=0.  # non zero default will be used
if [ $# -ge 5 ]; then
    norm=$5
fi
if [ $# -eq 6 ]; then
    xmin=$6
fi
if [ $func == "seltzer" ]; then
   E2=` echo $E2 | awk '{if($1 >= 10.) print 9.99; else print $1}' -`
fi
nbin=`awk  'END {nbin=int(log(E2/E1)/log(10.)/step)+1;print nbin}'   E1=$E1 E2=$E2 step=$step /dev/null`
Ek=$E1
if [ $func == "ps1" ]; then
    sel=1
else
    sel=3
fi
rm -f $work/fbrems.png
rm -f $work/fbrems.func

make -f DrawBremsFunc${sel}.mk
for f in $(seq 1 $nbin); do
  echo ${f}-th energy=$Ek 
  echo  $media $func $Ek $norm $xmin | ./drawbrems${sel}.out >> $work/fbrems.func
  echo " " >> $work/fbrems.func
  Ek=`awk 'END {Ek=Ek*10.**step;print Ek}' Ek=$Ek step=$step /dev/null`
done

echo "Ek=${E1}" > $work/gnuplot.com
echo "step=${step}"  >> $work/gnuplot.com
echo "media=\"${media}\""  >> $work/gnuplot.com
cat fbrems.gp >> $work/gnuplot.com
echo  set output  \"$work/fbrems.png\" >>$work/gnuplot.com
echo  rep  >>$work/gnuplot.com

gnuplot $work/gnuplot.com
