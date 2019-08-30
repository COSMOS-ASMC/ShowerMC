#!/bin/bash
cat <<EOF
a)  epos-lhc-v3700 (usable even A>56)
b)  epos-lhc-v3400 
c)  epos1.99
other) exit

Note:
 If projectile/target contains A>56,
 epos1.99/epos-lhc-v3400 cannot be used (say, in EPICS)            

 Which model ? Enter a,  b, c ..
EOF
read ans
if [ -z $ans ]; then
    exit
fi
source $COSMOSTOP/Scrpt/setarch.sh 

if [ $ans == "a" ]; then
    echo "epos-lhc-v3700 is being made to be selectable"
    $COSMOSTOP/Scrpt/ar-d-epos3400.sh
    $COSMOSTOP/Scrpt/ar-d-epos199.sh
   (cd $COSMOSTOP/Import/; rm -f EPOS; ln -s ./Hidden/EPOS/epos-lhc-v3700 EPOS;  cd $COSMOSTOP/Import/Hidden/EPOS/epos-lhc-v3700;    make veryclean; make) 

elif [ $ans == "b" ]; then
   echo "epos-lhc-v3400 is being made to be selectable"
    $COSMOSTOP/Scrpt/ar-d-epos3700.sh
    $COSMOSTOP/Scrpt/ar-d-epos199.sh
   ( cd $COSMOSTOP/Import/; rm -f EPOS; ln -s ./Hidden/EPOS/epos-lhc-v3400 EPOS;  cd EPOS;    make veryclean; make) 
elif [ $ans == "c" ]; then
    echo "epos199 is being made to be selectable"
    $COSMOSTOP/Scrpt/ar-d-epos3700.sh
    $COSMOSTOP/Scrpt/ar-d-epos3400.sh
    (cd $COSMOSTOP/Import/; rm -f EPOS; ln -s ./Hidden/EPOS/epos199  EPOS;  cd EPOS;    make veryclean; make) 
else
    exit
fi
