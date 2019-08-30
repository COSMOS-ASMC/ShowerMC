#!/bin/csh  -f
#  Usage: goto directory where you stored -t.hist data files 
#  and
#       $COSMOSTOP/UserHook/mkLDD/Util/check.csh
#  When a file is checked, program wait for your 'Enter'
#  If you don't find at least several lines for  
#  time bins, the file is invalid
#
echo "If you don't find several lines of time bins"
echo "The results are suspected invalid"
echo "If you specified 'no reduced time' for histogram"
echo "The first several time bins should be negative for inclied shower"
echo "If not, the result is invalid"

@ n=0
foreach f(*-t.hist)
  echo $f
  awk -f $COSMOSTOP/UserHook/mkLDD/Util/check.awk $f
  set dummy=$<
  @ n++
end
echo $n files checked

