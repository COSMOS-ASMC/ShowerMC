#!/bin/bash
# assume the input is raw srim output for stopping power.
#        dE/dx must be in MeV/(mg/cm2) = GeV/(g/cm2)
#        input file name is e.g, Zxx-SCIN.srim
#        where xx is the charge (atomic no.)
if [ $# != 3 ]; then
    cat <<EOF
    Usage: ./convSrim.sh Name forWhat ouputdir
       Name: media name such as SCIN
       forWhat:  1 --> for betagamma vs  dE/dx/Z**2
                 2 --> for Epics use  Ek vs dE/dx 
       outputDir;  output file is put there
EOF
    exit
fi

media=$1
forWhat=$2

for f in Z*-${media}.srim; do
    echo $f
    Z=`echo $f | awk '{i=index($0,"-"); Z=substr($0,2,i-2);printf("%02d\n", Z); }'`
    if [ $2 == 1 ]; then
	awk -f ./convSrim.awk Zin=$Z $f > $3/${media}-srim.$Z
    elif [ $2 == 2 ]; then
	awk -f ./convSrimForEpics.awk Zin=$Z $f > $3/srim.$Z
    fi
done
