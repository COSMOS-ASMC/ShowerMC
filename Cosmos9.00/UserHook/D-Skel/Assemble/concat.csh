#!/bin/csh -f
foreach f(*.dat-r)
    echo $f >> ended
    cat $f >> /Loft2/Work/kasahara/CosmosData/ForTA/p20/p20cos0.95/p20cos0.95.dat
    sleep 1
    rm $f
end

