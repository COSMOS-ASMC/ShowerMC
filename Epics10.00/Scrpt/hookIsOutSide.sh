#!/bin/bash
if [ -z "$EPICSTOP" ]; then
    echo EPICSTOP is not yet defined
    exit
fi

if [ ! -s ephook.f ]; then
    cat <<EOF
       ephook.f not exist in this directory=`pwd`
EOF
exit
fi
maininc=`awk '$1=="#include" && index($2,"ZepMain.f")>0 ' ephook.f`

if [ -f epicsfile ]; then
#   change "../../Data/Media" --> "$EPICSTOP/Data/Media"
   epicsfilemod=`awk '$1=="MediaDir" && index($2,"../Data/Media")>0 {print "yes"; exit}' epicsfile `
   if [ ! -z $epicsfilemod ]; then
       awk '$1=="MediaDir" && index($2,"../Data/Media")>0 {$1=" MediaDir"; $2="\"$EPICSTOP/Data/Media\""};{print}' epicsfile > temptemp.f
       mv temptemp.f epicsfile
   fi
fi


if [ ! -L $EPICSTOP/epics/ZepMain.f ]; then
    (cd $EPICSTOP/epics;  ln -s ../UserHook/main.f ZepMain.f)
    echo "link epics/ZepMain.f -->../UserHook/main.f" is made
fi


if  [ x"$maininc" != "x" ]; then
    echo "ok"
    exit
fi


if [ x"$maininc" == "x" ]; then
#     #include "../main.f" must be changed to #include "ZepMain.f"
   awk '$1=="#include" && index($2,"../main.f")>0 \
  { print "#include \"ZepMain.f\"" ;next};  \
    {print}'  ephook.f > temptemp.f
   mv temptemp.f ephook.f
   echo "ok"
   exit
fi
