#!/bin/csh -f
make clean
echo "select interaction model"
echo "1)  dpmjet3, qgsjet2, etc "
echo "2)  qgsjet1 "
echo "3)  sibyll  " 
set  ans=0
set  ans=$<

if ( $ans == 1 ) then
    echo "#define  OTHER" > Zintmodel.h
    make -f Gencol3.mk

else if ( $ans == 2 ) then
    if ( ! -f qgsjet01.f ) then
       if( -f $COSMOSTOP/Import/Test/QGS1/qgsjet01.f ) then
         ln -s $COSMOSTOP/Import/Test/QGS1/qgsjet01.f qgsjet01.f
       else
         echo "You need qgsjet program in " $COSMOSTOP/Import/Test/QGS1
	 echo "The name should be qgsjet01.f"
         echo "After installing, if the linker complains about "
	 echo "multiple definitions, then, "
         echo "rename them, say, adding _x to the names"
         echo "For the random number generator, it may be deleted."
	 echo "You need also SECTNU and QGSDAT01 there"
	 exit
       endif
    endif
    if ( ! -f qgs01init.f ) then
      ln -s $COSMOSTOP/Import/Test/QGS1/qgs01init.f qgs01init.f
    endif
    if ( ! -f SECTNU ) then
       ln -s $COSMOSTOP/Import/Test/QGS1/SECTNU SECTNU
    endif
    if ( ! -f QGSDAT01 ) then
       ln -s $COSMOSTOP/Import/Test/QGS1/QGSDAT01 QGSDAT01
    endif
    echo "#define  QGSJET1" > Zintmodel.h
    make -f Gencol3QGS1.mk
else if ( $ans == 3 ) then
    if ( ! -f sibyll2.1.f ) then
       if( -f $COSMOSTOP/Import/Test/Sibyll/sibyll2.1.f ) then
          ln -s $COSMOSTOP/Import/Test/Sibyll/sibyll2.1.f sibyll2.1.f
       else
         echo "You need sibyll in " $COSMOSTOP/Import/Test/Sibyll
	 echo "The name should be sibyl2.1.f"
	 exit
       endif
    endif
    if ( ! -f sibyllinit.f ) then
       ln -s $COSMOSTOP/Import/Test/Sibyll/sibyllinit.f sibyllinit.f
    endif
    echo "#define  SIBYLL" > Zintmodel.h
    make -f Gencol3Sibyll.mk
else
    echo invalid $ans
    exit
endif
echo " "
echo "************************************************************** "
echo "NOTE: You have to specify a corresponding model by IntModel in param"
echo " "
if ( $ans == 2 ) then
   echo 'IntModel cannot be "sibyll"'
   echo '  but can be "qgsjet1" "qgsjet2" "dpmjet3" etc'
else if ($ans == 3 ) then
   echo 'IntModel cannot be "qgsjet1"'
   echo '  but can be "sibyll" "qgsjet2" "dpmjet3" etc'
else
    echo 'IntModel cannot be "qgsjet1" nor "sibyll"'
    echo ' but can be "qgsjet2" "dpmjet3" etc'
endif
echo " "
