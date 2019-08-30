#!/bin/bash
if [ -z "$COSMOSTOP" ]; then
    echo COSMOSTOP is not yet defined
    exit
fi

if [ ! -s chook.f ]; then
    cat <<EOF
       chook.f not exist in this directory=`pwd`
EOF
exit
fi

maininc=`awk '$1=="#include" && index($2,"../cmain.f")>0 {print "ng";exit} ' chook.f`
echo "maininc=" $maininc

if [ ! -L $COSMOSTOP/cosmos/cmain.f ]; then
echo "make cmain link"
    (cd $COSMOSTOP/cosmos;  ln -s ../UserHook/cmain.f cmain.f)
    echo "link cosmos/cmain.f -->../UserHook/cmain.f" is made
fi

sed   "s/\.\/Zabsorb.h/Zabsorb.h/g"   $COSMOSTOP/UserHook/chookEabsorb.f > temp
mv temp  $COSMOSTOP/UserHook/chookEabsorb.f

if [ ! -L $COSMOSTOP/cosmos/chookEabsorb.f ]; then
    (cd  $COSMOSTOP/cosmos/; ln -s ../UserHook/chookEabsorb.f chookEabsorb.f )
fi
    

if [ ! -L $COSMOSTOP/cosmos/Zabsorb.h ]; then
   echo "make Zabsorb.h link"
    (cd $COSMOSTOP/cosmos;  ln -s ../UserHook/Zabsorb.h Zabsorb.h )
    echo "link cosmos/Zabsorb.h-->../UserHook/Zabsorb.h" is made
fi

if [ x"$maininc" == "xng" ]; then
#     #include "../cmain.f" must be changed to #include "cmain.f"
   awk '$1=="#include" && index($2,"../cmain.f")>0 \
  { print "#include \"cmain.f\"" ;next};  \
    {print}'  chook.f > temptemp.f
   mv temptemp.f chook.f
   echo "ok 2"
fi
Cereninc=`awk '$1=="#include" && index($2,"../ctemplCeren.f")>0 {print "ng";exit} ' chook.f`
if  [ x"$Cereninc" != "x" ]; then
    echo "ok 3"
fi
if [ x"$Cereninc" == "xng" ]; then
#     #include "../ctemplCeren.f" must be changed to #include "ctemplCeren.f"
   awk '$1=="#include" && index($2,"../ctemplCeren.f")>0 \
  { print "#include \"ctemplCeren.f\"" ;next};  \
    {print}'  chook.f > temptemp.f
   mv temptemp.f chook.f
   echo "ok 4"
   exit
fi

