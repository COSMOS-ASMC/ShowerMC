#!/bin/csh
#  From 'rvsrho.data' data created  (r vs 2pirrho)
#  by procLat
#  which processes  -r.hist data,
#  this makes r vs Rho for series of  fai, code.
#  and make a fitting of rho*2pi=a/r**(b+csqrt(r)) 
#  and put coefficients a,b,c  in ./fitted.data.
#  ./fitted.data must not exist beforehand
if ( $#argv < 2 || $#argv > 4 ) then
    cat <<EOF
    Usage: FitRvsRho.csh  inputrvsrhofile  layer  {age cogd}
    inputrvsrhofile:  path to the rvsrho.data file
    (made by procLat.f) 
    layer: input file is for 'layer' in  .hyb file"
       layer is not used at present.
    age: optional age
    cogd: optional cog depth   
EOF
    exit 1
endif

if ( -f ./fitted.data ) then
    echo "./fitted.data  exist"
    exit 1
endif

set rrhofile=$1
set layer=$2
#   ********************
#    ptcl code
set codeA=(1 2 3 4)
#     code symbol
set codeN=(g e m h)
#     fai indexes
set faiA=(1 2 3 4 5 6 7 8 9 10 11 12)
#     actual value in deg.
set faiV=(0 30 60 90 120 150 180 210 240 270 300 330)
#     *******************
set ARCH=`../setarch.sh`
#   get fitting executable from baseInfo
set fitter=`awk '$1=="fitter" {print $2}' ./baseInfo`${ARCH}

echo "fitter is $fitter"

#  next could be simply ./
set head=`pwd`
if( $#argv == 4 ) then
    echo $3 $4 > ./fitted.data
endif

@ codec=0
foreach code($codeA)
    @ codec++
    @ faic=0 
    foreach fai($faiA)
	@ faic++
#                   layer
        echo "l " $layer  $code $fai  >> ./fitted.data
        awk -f $head/getLatAtFai.awk fai=$fai code=$code $rrhofile  |  $fitter  $code c  >> ./fitted.data
       if ( $status != 0 ) then
          echo "fitting failed"
          exit 1
       endif
    end
 end
echo "l 0 0 0" >> ./fitted.data
