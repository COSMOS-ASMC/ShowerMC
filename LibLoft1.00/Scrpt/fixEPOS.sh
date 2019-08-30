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
source $LIBLOFT/Scrpt/setarch.sh 
# delete all possible epos related object in libloft.a
echo "---Don't worry about many messages like:"
echo "   ar: cepos.o: not found in archive--- "
$LIBLOFT/Scrpt/ar-d-epos.sh


if [ $ans == "a" ]; then
    echo "epos-lhc-v3700 is being made to be selectable"
#    $LIBLOFT/Scrpt/ar-d-epos3400.sh
#    $LIBLOFT/Scrpt/ar-d-epos199.sh
   (cd $LIBLOFT/Had/Import/; rm -f EPOS; ln -s ./Hidden/EPOS/epos-lhc-v3700 EPOS;  cd $LIBLOFT/Had/Import/Hidden/EPOS/epos-lhc-v3700;    make veryclean; make) 

elif [ $ans == "b" ]; then
   echo "epos-lhc-v3400 is being made to be selectable"
#   $LIBLOFT/Scrpt/ar-d-epos3700.sh
#    $LIBLOFT/Scrpt/ar-d-epos199.sh
   ( cd $LIBLOFT/Had/Import/; rm -f EPOS; ln -s ./Hidden/EPOS/epos-lhc-v3400 EPOS;  cd EPOS;    make veryclean; make) 
elif [ $ans == "c" ]; then
    echo "epos199 is being made to be selectable"
#    $LIBLOFT/Scrpt/ar-d-epos3700.sh
#    $LIBLOFT/Scrpt/ar-d-epos3400.sh
    (cd $LIBLOFT/Had/Import/; rm -f EPOS; ln -s ./Hidden/EPOS/epos199  EPOS;  cd EPOS;    make veryclean; make) 
else
    exit
fi
