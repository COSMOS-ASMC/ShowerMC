#!/bin/csh
foreach f(`awk '{print $1}' fini`)
   echo $f
   foreach g(`grep $f jobs`)
     echo $g
     rsync -e ssh -avz ${f}:/tmp/kasahara2/"*"$g"*" /tmp/kasahara2/ >& error
   end
 end
