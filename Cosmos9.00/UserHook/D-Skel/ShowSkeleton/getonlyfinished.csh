#!/bin/csh 
#  adjust substr positions below;
#  ThinHost  must be adjusted  by seeing hostnum
# use reduceEach.bin.atsite.csh and later  use getonlyReduced.csh
#
rm -f finished
foreach f( *.err )
  set x=`tail -4  $f | grep comp`
  if( "x$x" != "x" ) then
    echo $f >> finished
  endif
end
#
awk '{print substr($0,33,8), substr($0,42,3)}' finished > hostnum



foreach hn ( "`cat hostnum`" )
 set h=`echo $hn | awk '{print $1}' `
 set n=`echo $hn | awk '{print $2}' `
 rsync -e ssh -avz --exclude "*".dat  ${h}:/tmp/kasahara2/"*"$n"*" /tmp/kasahara2/
end

