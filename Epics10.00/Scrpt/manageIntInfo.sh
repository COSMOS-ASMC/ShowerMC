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

modified=`awk '$1=="#include" && index($2,"epUI.f")>0 ' ephook.f`


if [ ! -s epUI.f ] ; then
    cp $EPICSTOP/UserHook/epUI.f ./
    echo "epUI.f has been copied to this directory"
fi

if [ x"$modified" != "x" ]; then
    echo "ok"
    exit
fi

if  [ x"$modified" == "x" ]; then
#   find  #include "../main.f" 
   awk '$1=="#include" && index($2,"../main.f")>0 \
    { print ; system("cat $EPICSTOP/UserHook/INTINFOinclude"); next };  \
    {print}'  ephook.f > temptemp.f
   cmp -s ephook.f temptemp.f
   if [ $? -ne 0 ]; then
       mv temptemp.f ephook.f
       echo "ok-2"
       exit
   fi
   rm -f temptemp.f 
fi

if [ x"$modified" == "x" ]; then
#     find  #include ZepMain.f"
   awk '$1=="#include" && index($2,"ZepMain.f")>0 \
  { print ;system("cat $EPICSTOP/UserHook/INTINFOinclude"); next };  \
    {print}'  ephook.f > temptemp.f
   mv temptemp.f ephook.f
   echo "ok"
fi
