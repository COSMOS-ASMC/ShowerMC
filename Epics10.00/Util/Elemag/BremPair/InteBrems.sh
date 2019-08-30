#!/bin/bash
work=Work
if [ $#  -lt  8 ] || [ $# -gt 10 ]; then
cat <<EOF
  This integrates  a specified brems function in a given
  energy region even if it is outsied of the proper
  region. 

 Usage: ./InteBrems.sh  media func how E1 E2 step outfile  [norm
           LPMeffect Flpm]
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
   outfile: output file name

    norm: normalization.   (default =1)
           1--> / r.l 
           2--> mb/ingredient
           3--> /(g/cm^2)
           4--> /cm
           5--> area normlization.  Integral (xmin,xmax) of dsigma/dx=1
   LPMeffect: (default is t). If f, LPM is neglected
   Flpm:    (default is 1). (must be >=1).  LPM works when Ee
        > Flpm*(minimum LPM energy) and LPMeffect=t
EOF
exit
fi
media=$1;  
func=$2;
how=$3
echo "media=" $media " func=" $func "how=" $how
E1=$4; E2=$5; step=$6;
echo "E1=" $E1 " E2=" $E2 " step=" $step

outfile=$7 

if [ $# -ge 8 ]; then
    norm=$8
else
    norm=1 
fi


if [ $# -ge 9 ]; then
    LPMeffect=$9
else
    LPMeffect="t"
fi

if [ $# -ge 10 ]; then
    Flpm=$10
else
    Flpm=1
fi

echo "norm=" $norm "LPMeffect=" $LPMeffect "Flpm=" $Flpm  
echo " "

echo "sleeing 3 seconds"
sleep 3


Ek=$E1
if [ $func == "ps3" ]; then
    sel=3
else
    sel=1
fi

make clean
make -f InteBrems${sel}.mk
if [ $? != 0 ]; then
    echo compile error
    exit 1
fi

echo  $media $func  $how  $E1 $E2 $step $norm $LPMeffect $Flpm | ./Intebrems${sel}.out  > $outfile

