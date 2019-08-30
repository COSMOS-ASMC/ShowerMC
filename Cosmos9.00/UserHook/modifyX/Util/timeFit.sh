#!/bin/bash
#  This is a dedicated routine for LDD.
#  This is to get coefficients a,b,c in T10=a*r**(b+c*log(r));
#  by calling  procTime and  mkRvsTandFit.csh
#  and store them in a given directory.

#  next is not used. mu0 is obtained from 
#  mu0=`awk '$1=="mu0" '{print $2}' ./baseInfo`  

if [  $# -ne 4  ] ; then
    cat << EOF 
     Usage: ./timeFit.sh timehist cosz layer# mu
   
  timehist: -t.hist file path 
      cosz:  cos of the  primary zenith angle
    layer#:  a layer number in .hyb data where -t.hist data
             was consturcted
        mu:  Molire unit length (m) at the layer
EOF
echo "# of arguments is " $#
exit 1    
fi
ARCH=`./setarch.sh`
tfile=$1
cosz=$2
ly=$3
mu=$4
smooth=`awk '$1=="smooth" {print $2}' ./baseInfo`
if [ x$smooth = "x" ]; then
    echo "no smooth def. in baseInfo"
    exit 1
fi
reduced=`awk '$1=="reducedT" {if($2=="yes") print 1; else print 0 }' ./baseInfo`
if [ x$reduced = "x" ]; then
    echo "no reducedT def. in baseInfo"
    exit 2
fi
# judge input if ascii hist or binary hist
ascii=(`file -b $tfile`) 
if [ ${ascii[0]} = "ASCII" ] ; then
    format=1
else
    format=2
fi
#  next 2 is to get upto T10%
#  next 0 is for hist  from mkLDD (contains hist for core region)
 ./procTime${ARCH}  $format $reduced 2  0 $smooth  $tfile > ./tf.data

if [ $? -ne 0 ]; then
     echo "procTime${ARCH} failed"
     exit 1
else
    echo  "procTime${ARCH} ended" 
fi

 ./mkRvsTandFit.csh $mu  $mu $cosz tf.data $ly
 if [ $? -ne 0 ]; then
     echo mkRvsTandFit failed
     exit 1
 fi 

  echo  "mkRvsTandFit ended"
