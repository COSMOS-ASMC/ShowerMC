#!/bin/csh
if( $#argv != 4 ) then
    cat <<EOF
    ./selectDataForplot.csh  faiIdx minAge maxAge inputFile 
	InputFle is the output from getRhoAtFewR.csh, i.e, 
    rhoAtRs.data
EOF
    exit 1
endif
set codeA=(gamma elec muon hadron)
@ n=1
while ($n < 5) 
#    echo $n
    awk -f selectDataForPlot.awk code=$n fai=$1 minAge=$2 maxAge=$3 $4 >  $codeA[$n]S${2}-$3.data
   ./mkHist.csh $codeA[$n]S${2}-$3.data
    @ n++
end


