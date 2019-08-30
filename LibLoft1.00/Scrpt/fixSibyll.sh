#!/bin/bash
cat <<EOF
a)  sibyll2.3c
b)  sibyll2.1
other) exit

Note:
 If projectile contains A>56,
 sibyll cannot be used (say, in EPICS)            
 It can be used only for proton and "Air" target.
   (Ar is not contained)
   It is dangerouse to use it even for N2 or O2 target
 Which model ? Enter a,  b  or ..
EOF
read ans
if [ -z $ans ]; then
    exit
fi
source $LIBLOFT/Scrpt/setarch.sh
# remove all sibyll related object files in libloft.a
 echo "---Don't worry about many messages like: " 
 echo "     ar: csibyll.o: not found in archive---"  
 $LIBLOFT/Scrpt/ar-d-sibyll.sh

if [ $ans == "a" ]; then
    echo "sibyll2.3c is being made to be selectable"
#    $LIBLOFT/Scrpt/ar-d-sibyll2.1.sh
   (cd $LIBLOFT/Had/Import/; rm -f Sibyll; ln -s ./Hidden/Sibyll/sibyll2.3c Sibyll;  cd $LIBLOFT/Had/Import/Hidden/Sibyll/sibyll2.3c;    make veryclean; make) 

elif [ $ans == "b" ]; then
   echo "sibyll2.1 is being made to be selectable"
#    $LIBLOFT/Scrpt/ar-d-sibyll2.3c.sh
   ( cd $LIBLOFT/Had/Import/; rm -f Sibyll; ln -s ./Hidden/Sibyll/sibyll2.1 Sibyll;  cd Sibyll;    make veryclean; make) 
else
    exit
fi
