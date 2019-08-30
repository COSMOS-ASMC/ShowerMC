#!/bin/csh -f
set xxx=`ps aux | grep kasahara | awk '$0 !~ "grep"' | awk '{print $2;exit}' `
echo "xxx=" $xxx
if( x$xxx == "x" ) exit
echo job $xxx found
kill -9 $xxx

