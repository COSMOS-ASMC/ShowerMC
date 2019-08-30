#!/bin/bash
#  This is a dedicated routine for LDD.
#  This is to get coefficients a,b,c in 
#  rho*2pir=a/r**(b+c*sqrt(r));
#  by calling  procLat and  latFit program
#  and store them in a given directory.


if [ $# -lt 2  -o  $# -gt 4   ] ; then
    cat << EOF 
     Usage: ./latFit.sh rhist  layer# {age cogd}
   
  rhist: -r.hist file path 
    layer#:  a layer number in .hyb data where -r.hist data
             was consturcted
     age:   optional age
     cogd:  optional depth
EOF
echo "# of arguments is " $#
exit 1    
fi
ARCH=`../setarch.sh`
rfile=$1
ly=$2
age=$3
cogd=$4
# judge input if ascii hist or binary hist
ascii=(`file -b $rfile`) 
if [ ${ascii[0]} = "ASCII" ] ; then
    format=1
else
    format=2
fi
 ./procLat${ARCH}  $format  $rfile > ./rvsrho.data

if [ $? -ne 0 ]; then
     echo "procLa${ARCH} failed"
     exit 1
else
    echo  "procLat${ARCH} ended" 
fi

 ./FitRvsRho.csh rvsrho.data $ly $age $cogd
 if [ $? -ne 0 ]; then
     echo FitRvsRho.csh failed
     exit 1
 fi 

  echo  "FitRvsRho.csh ended"
