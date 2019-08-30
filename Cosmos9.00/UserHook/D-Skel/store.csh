#!/bin/tcsh
@ n=0
foreach f(`awk '{print $1}' ThinHosts `)
   echo $f
   cat /tmp/kasahara2/*$f.dat-r >> /Loft2/Work/kasahara/CosmosData/ForTA/p20/dpm/p20cos0.975eqsd/p20cos0.975E900-45.dat
    @ n++
    if( $n == 45 ) exit
end
