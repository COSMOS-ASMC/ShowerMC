#!/bin/csh
#  This makes r vs T for series of  fai, code and percentage
#  
set head=$COSMOSTOP/UserHook/mkLDD/Util

if( $#argv != 5 ) then
    cat <<EOF
    Usage: mkRvsT.csh  m.u0  m.u cosz  inputTFdata outputfile
    m.u0: Moliere unit length (in m) of the  assmued
          observation site: 86 for 875 g/cm2
    m.u:  is the Moliere unit length (in m) of the height
          where the input data is constructed. m.u could be m.u0
    cosz: cosine of the 1ry zenith angle. If reduced time output is
          not wanted (or already in reduced time), 
          **** use a vaule <=0. independent of actual cos****
          Note: cosz=1-->reduced time and actual time coinside
	    However, it is not suited for curve fitting and
	    we force to use "reduced time" using cosz=0.999. 
          Normally for cos=1.0, put 1.0   
  inputTFdata: such as tf.data made by procTime.f
  outputfile: dir/basic file name of output. 
         which will contain r(m.u) and time(ns)
         time in the output data is reduced time (if cosz> 0).
         r is scaled to the depth of which Moliere unit is m.u0
         (r in the file is converted to r*mu/mu0 and used as r
	 in the (r, time).
           output file name will be automatically created like 
         (suppose it is /tmp/LDD)
         /tmp/LDDgF1p50.data etc (for gamma, faiindex = 1, percentage=50%)
  ** NOTE:  the script uses $head/mkRvsTcond.csh for what fai angles,
            particle codes, percentage of tf data are to be used.
EOF
    exit
endif

source $head/mkRvsTcond.csh


@ codec=0
foreach code($codeA)
    @ codec++
    @ faic=0 
    foreach fai($faiA)
	@ faic++
	@ pcc=0
	foreach pc($percentA)
	    @ pcc++
#                  nth=..  mu=.. mu0=.. fai=.. cosz=..
	    awk -f $head/getTimeAtFai.awk  fai=$fai code=$code $4  |  awk -f $head/getTimeAtFrac.awk nth=$pc mu0=$1 mu=$2 fai=$faiV[$faic] cosz=$3 > $5$codeN[$codec]F${fai}p$percentV[$pcc].data
        end
     end
end

