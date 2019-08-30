#!/bin/csh
#  From 'tf' data created by procTime which processes  -t.hist data,
#  this makes r vs T for series of  fai, code and percentage
#  and make a fitting of T=a*r**(b+clog(r)) to (r,T) 
#  and put coefficients a,b,c  in ./fitted.data.
#  ./fitted.data must not exist beforehand
if ( $#argv != 5 ) then
    cat <<EOF
    Usage: mkRvsTandFit.csh  mu0  mu cosz  inputTFfile  layer

     m.u0: Moliere unit length (in m) of the observation site
           86 or 87 for 875 g/cm2
      m.u: Moliere unit length (in m) of the height where the input
           'tf' data is constructed. For LDD,  mu should be mu0 
     cosz: cosine of the 1ry zenith angle. 
    inputTFfile:  path to the tf file (such as tf.data made by procTime.f) 
    layer: input file is for 'layer' in  .hyb file"
EOF
    exit 1
endif

if ( -f ./fitted.data ) then
    echo "./fitted.data  exist"
    exit 1
endif

set mu0=$1
set mu=$2
set cosz=$3
set tffile=$4
set layer=$5
#  if -t.hist is reduced time, make cosz negative 
#  however, for cosz=1 (reduce or non-reduce is the same)
#  we force to use reduced time using cosz=0.999
set reducedtime=`awk '$1=="reducedT" {print $2}' ./baseInfo`

echo "reducedtime=$reducedtime"

if ( $reducedtime == "yes" ) then
    set needredt=`echo $cosz | awk '{ if( $1==1.0) print "yes"; else print  "no"}'`
    
#    echo "needredt= $needredt"

    if ( $needredt == "no" ) then
	set cosz=-1.
    endif
else
#            not used 
  set  needredt="yes"
endif

# echo "reducedtime=" $reducedtime " needredt="$needredt " cosz="$cosz
#  cosz<0 means we need not convert time to reduced time.

#   ********************
#    ptcl code
set codeA=(1 2 3 4)
#     code symbol
set codeN=(g e m h)
#     fai indexes
set faiA=(1 2 3 4 5 6 7 8 9 10 11 12)
#     actual value in deg.
set faiV=(0 30 60 90 120 150 180 210 240 270 300 330)
#     T10 % only; 3rd column 
set percentA=(3) 
#       is 10 %
set percentV=(10)
#     *******************
set ARCH=`./setarch.sh`
#   get fitting executable from baseInfo
set fitter=`awk '$1=="fitter" {print $2}' ./baseInfo`${ARCH}

echo "fitter is $fitter"

#  next could be simply ./
set head=`pwd`

@ codec=0
foreach code($codeA)
    @ codec++
    @ faic=0 
    foreach fai($faiA)
	@ faic++
	@ pcc=0
	foreach pc($percentA)
	    @ pcc++
#                     layer
	    echo "l " $layer  $code $fai >> ./fitted.data
	    awk -f $head/getTimeAtFai.awk  fai=$fai code=$code $tffile  |  awk -f $head/getTimeAtFrac.awk nth=$pc mu0=$mu0 mu=$mu fai=$faiV[$faic] cosz=$cosz |  $fitter  $code c  >> ./fitted.data
	    if ( $status != 0 ) then
	       echo "fitting failed"
	       exit 1
 	    endif
        end
     end
end
echo "l 0 0 0" >> ./fitted.data
