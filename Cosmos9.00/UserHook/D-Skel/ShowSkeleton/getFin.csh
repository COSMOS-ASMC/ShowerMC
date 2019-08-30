#!/bin/csh
awk '{printf("%s.%s\n", $1,$2)}' hostnum > jobs
foreach f(`awk '{print $1}' hostnum`)
   echo $f
   foreach g(`grep $f jobs`)
     echo $g
     rsync -e ssh -avz ${f}:/tmp/kasahara/"*"$g"*" /tmp/kasahara/ >& error
   end
 end
