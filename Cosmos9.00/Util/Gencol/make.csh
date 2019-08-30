#!/bin/csh -f
make clean
cat <<EOF
select interaction model
1)  dpmjet3, qgsjet2, epos, sibyll etc 
2)  qgsjet1 (not useable now)
Note: 
qgsjet2:
  To switch qgsjetII-04<-->qgsjetII-03, use fixIntModel.sh
  at any place. Default at release time is 04. 
  Note: 03 needs qgsdat-II-03.ascii in Cosmos/Data/QGS 
        04 needs qgsdat-II-04 in
                     Cosmos/Import/Hidden/QGS/qgsjetII-04
        They are too large to include in each Cosmos version.
       (to get them, visit http://cosmos.n.kanagawa-u.ac.jp/)
epos:
 To switch 3 versions of EPOS, use fixIntModel.sh
 To see which version is currently usable, use showIntModel.sh .

 qgsjet1:
    At present, this does not work; (name collision with EPOS) etc
    Use older versions (say, Epics9.08 + Cosmos7.622; however,
    qgsjet1 with Diffraction information cannot be used (reason unknown;
    there should be no difference in the produced particle spectra).
EOF
set  ans=0
set  ans=$<
echo your ans is $ans
if ( x$ans == "x" ) then
    exit
endif

cat<<EOF
Don't worry next, if you don't know artificial modifiation.
  Do you try to artificial modification of pi/K/eta X?
  If yes, program in Cosmos/UserHook/modifyX/csoftenPiK.f will be used
  Enter y, if so
EOF
set  ans2=0
set  ans2=$<
if ( x$ans2 == "xy" ) then
    echo "#define MODIFYX" > Zintmodel.h
    
    if ( ! -f csoftenPiK.f ) then
       echo "csoftenPiK.f is not seen here "
       if ( -l csoftenPiK.f ) then
           echo "csoftenPiK.f is link; so I'm rm it and make link"
  	   rm -f csoftenPiK.f
       endif
       set x="need"
    else
#       even if csoften.. exists, it may be a link to old Cosmos; must renew
       set x=`ls -l csoftenPiK.f | grep  $COSMOSTOP | awk 'NF>0 {print $1}'`
    endif
    if ( x$x == "x" || $x == "need" ) then
       rm -f csoftenPiK.f
       if ( -f $COSMOSTOP/UserHook/modifyX/csoftenPiK.f ) then
         echo "make a link to csoftenPiK.f"
         ln -s $COSMOSTOP/UserHook/modifyX/csoftenPiK.f  csoftenPiK.f
       else
         echo "$COSMOSTOP/UserHook/modifyX/csoftenPiK.f non existent"  
         exit
       endif
    endif

    if ( ! -f softenparam.dat ) then
       echo "softenparam.dat is not seen here "
       if ( -l softenparam.dat ) then
         echo "softenparam.dat is link; so I'm removing it and make link"
         rm -f softenparam.dat
       endif
       set x="need"
    else
       set x=`ls -l csoftenparam.dat | grep  $COSMOSTOP | awk 'NF>0 {print $1}'`
    endif
    if( x$x == "x" || $x == "need" ) then
       rm -f softenparam.dat 
       if ( -f $COSMOSTOP/UserHook/modifyX/softenparam.dat ) then
          echo "making a link to softenparam.dat"
          ln -s $COSMOSTOP/UserHook/modifyX/softenparam.dat softenparam.dat 
       else
           echo "$COSMOSTOP/UserHook/modifyX/softenparam.dat non existent"  
	   exit
       endif
    endif
else
    rm -f Zintmodel.h
    touch Zintmodel.h
endif

if ( $ans == 1 ) then
    echo "#define  OTHER" >> Zintmodel.h
    make -f Gencol.mk

else if ( $ans == 2 ) then
    if ( ! -f qgsjet01.f ) then
       echo "qgsjet01.f is not seen"
       if ( -l qgsjet01.f ) then
           echo "qgsjet01.f is link; so rm it and make link"
  	   rm qgsjet01.f
       endif
       echo COSMOSTOP is $COSMOSTOP
       if(  -f $COSMOSTOP/Import/Test/QGS1/withDiffFlag/qgsjet01.f ) then
         echo qqgjet01.f is found in $COSMOSTOP/Import/Test/QGS1/withDiffFlag/
         ln -s $COSMOSTOP/Import/Test/QGS1/withDiffFlag/qgsjet01.f qgsjet01.f
       else
         echo qqgjet01.f is NOT found in $COSMOSTOP/Import/Test/QGS1/withDiffFlag
         echo "After installing, if the linker complains about "
	 echo "multiple definitions, then, "
         echo "rename them, say, adding _x to the names"
         echo "For the random number generator, it may be deleted."
	 echo "You need also SECTNU and QGSDAT01 there"
	 exit
       endif
    endif
    echo "compiling qgsjet1"
    if ( ! -f qgs01init.f ) then
      rm -f qgs01init.f 
      ln -s $COSMOSTOP/Import/Test/QGS1/qgs01init.f qgs01init.f
    endif
    if ( ! -f SECTNU ) then
       rm -f SECTNU
       ln -s $COSMOSTOP/Import/Test/QGS1/SECTNU SECTNU
    endif
    if ( ! -f QGSDAT01 ) then
       rm -f  QGSDAT01 
       ln -s $COSMOSTOP/Import/Test/QGS1/QGSDAT01 QGSDAT01
    endif
    echo "#define  QGSJET1" >> Zintmodel.h
    make -f GencolQGS1.mk
else
    echo invalid $ans
    exit
endif
echo " "
echo "************************************************************** "
echo "NOTE: You have to specify a corresponding model by IntModel in param"
echo " "
cat <<EOF

If you are going to use dpmjet3, the energy in dpmjet.inp
must be little bit larger than E =( s-m**2)/(2m)
where s is (P1+P2)**2 =(E1+E2)**2- (\vec p1+ \vec p2)**2
where \vec p1 and \vec p2 are those given in param file
  If you are not sure about E value, run the program, 
  then, you will find it
EOF
