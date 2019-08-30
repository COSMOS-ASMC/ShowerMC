#!/bin/bash
work=Work
if [ $#  -ne  6 ] ; then
cat <<EOF
 Usage: $0  media xmin E1 E2 step  norm
   This computes total Brems Xsection in the 
   Seltzer-Berger numerical cross-section region.
   Integration is performed from xmin to xmax

   media: such as Fe BGO Air ...
   xmin = Egmin/Ee   >0.  xmax= 1 (not Egmin/(Ee-Me))
          To see the consistency with the cross-section
          at the Screening region, we may set xmin=10^-3 or so;
          at lower xmin, the difference is not so important.
   E1 E2: kinetic energy of electron; form E1 to E2
          E1 > 10^-6 GeV E1 <E2 < 10GeV
    step:  log10 step of energy
    norm:  normalization.   
           1--> / r.l 
           2--> mb/ingredient
           3--> /(g/cm^2)
           4--> /cm
           5--> area normlization.  Integral (xmin,1) of dsigma/dk=1
EOF
exit
fi
   
media=$1; xmin=$2;  E1=$3; E2=$4; step=$5; norm=$6;



rm -f $work/bremXS
make -f getTotalXSofSeltzerBremsFunc.mk

 echo  $media $xmin $E1 $E2 $step  $norm  | ./getTotalXSofSeltzerBremsFunc.out >> $work/bremXS


echo "output is in" $work/bremXS
