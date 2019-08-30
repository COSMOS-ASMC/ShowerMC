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
source $LIBLOFT/Scrpt/setarch.sh 
  (cd $LIBLOFT; \
      ar -d ./lib/$ARCH/libloft.a qgsjet-II-03.o; \
      ar -d ./lib/$ARCH/libloft.a cQGSjet.o; \
      ar -d ./lib/$ARCH/libloft.a qgsjet-II-04.o; 
      )

      if [ $ans == "a" ]; then
    echo "qgsjetII-04 is being made to be selectable"
  (cd $LIBLOFT; \
#      ar -d ./lib/$ARCH/libloft.a qgsjet-II-03.o; \
#      ar -d ./lib/$ARCH/libloft.a cQGSjet.o; \
      cd ./Had/Import/; rm -f QGS; ln -s ./Hidden/QGS/qgsjetII-04 QGS;  cd QGS;    make veryclean; make) 
###   echo "For qgsjet1 in Gencol, psran is needed so GencolQGS1.mk there is changed"
###   cp $LIBLOFT/Util/Gencol/GencolQGS1.mkWithQGSII-04  $LIBLOFT/Util/Gencol/GencolQGS1.mk

elif [ $ans == "b" ]; then
   echo "qgsjetII-03 is being made to be selectable"
  (cd $LIBLOFT; \
#      ar -d ./lib/$ARCH/libloft.a qgsjet-II-04.o; \
#      ar -d ./lib/$ARCH/libloft.a cQGSjet.o; \
      cd ./Had/Import/; rm -f QGS; ln -s ./Hidden/QGS/qgsjetII-03 QGS;  cd QGS;    make veryclean; make) 

###   echo "For qgsjet1 in Gencol, psran is not  needed so GencolQGS1.mk there is changed"
###   cp $LIBLOFT/Util/Gencol/GencolQGS1.mkWithQGSII-03  $LIBLOFT/Util/Gencol/GencolQGS1.mk
else
    echo your input, $ans, invalid
    exit
fi
