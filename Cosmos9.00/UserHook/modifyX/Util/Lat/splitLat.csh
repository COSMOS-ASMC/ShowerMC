#!/bin/csh
#  This makes r vs dn/dr  for series of  fai, code.
#  
set head=$COSMOSTOP/UserHook/mkLDD/Util/Lat

if( $#argv != 2 ) then
    cat <<EOF
    Usage: splitLat.csh   inputLatdata outputfile
  inputLatdata: such as Lat.data made by procLat.f
  outputfile: dir/basic file name of output. 
         which will contain r(m.u) and dN/dr
         (suppose it is Work/LDDLat)
         Work/LDDLatgF1.data etc (for gamma, faiindex = 1)
  ** NOTE:  the script uses $head/splitLatcond.csh for what fai angles,
            particle codes of Lat.data are to be used.
EOF
    exit
endif

source $head/splitLatcond.csh


@ codec=0
foreach code($codeA)
    @ codec++
    @ faic=0 
    foreach fai($faiA)
	@ faic++
#                  nth=..  mu=.. mu0=.. fai=.. cosz=..
       awk -f $head/getLatAtFai.awk  fai=$fai code=$code $1 > $2$codeN[$codec]F${fai}.data
     end
end
