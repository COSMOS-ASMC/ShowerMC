#!/bin/csh 
# ThinHost is assumed




foreach hn ( "`cat ../ThinHost`" )
 set n=`echo $hn | awk '{print $1}' `
 set h=`echo $hn | awk '{print $2}' `
# rsync -e ssh -avz  ${h}:/tmp/kasahara2/"*"$n"*-r" /tmp/kasahara2/
 rsync -e ssh -avz  ${h}:/tmp/kasahara2/"*"$n"*-r" /Work1/temp/
end

