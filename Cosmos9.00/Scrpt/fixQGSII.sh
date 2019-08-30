#!/bin/bash

cat <<EOF
a)  qgsjetII-04  (LHC tuned)
b)  qgsjetII-03
Which model ? Enter a or b.
EOF
read ans
if [ -z $ans ]; then
    exit
fi
source $COSMOSTOP/Scrpt/setarch.sh 

if [ $ans == "a" ]; then
    echo "qgsjetII-04 is being made to be selectable"
  (cd $COSMOSTOP; \
      ar -d ./lib/$ARCH/libcosmos.a qgsjet-II-03.o; \
      ar -d ./lib/$ARCH/libcosmos.a cQGSjet.o; \
      cd ./Import/; rm -f QGS; ln -s ./Hidden/QGS/qgsjetII-04 QGS;  cd QGS;    make veryclean; make) 
###   echo "For qgsjet1 in Gencol, psran is needed so GencolQGS1.mk there is changed"
###   cp $COSMOSTOP/Util/Gencol/GencolQGS1.mkWithQGSII-04  $COSMOSTOP/Util/Gencol/GencolQGS1.mk

elif [ $ans == "b" ]; then
   echo "qgsjetII-03 is being made to be selectable"
  (cd $COSMOSTOP; \
      ar -d ./lib/$ARCH/libcosmos.a qgsjet-II-04.o; \
      ar -d ./lib/$ARCH/libcosmos.a cQGSjet.o; \
      cd ./Import/; rm -f QGS; ln -s ./Hidden/QGS/qgsjetII-03 QGS;  cd QGS;    make veryclean; make) 

###   echo "For qgsjet1 in Gencol, psran is not  needed so GencolQGS1.mk there is changed"
###   cp $COSMOSTOP/Util/Gencol/GencolQGS1.mkWithQGSII-03  $COSMOSTOP/Util/Gencol/GencolQGS1.mk
else
    echo your input, $ans, invalid
    exit
fi
