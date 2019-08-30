#!/bin/csh
foreach f(`awk '{print $2}' allHosts`)
    echo $f
    rsync -e ssh -av {$f}:/tmp/kasahara/"*"dat /Loft2/Work/kasahara/CosmosData/ForTA/ 
end



  
