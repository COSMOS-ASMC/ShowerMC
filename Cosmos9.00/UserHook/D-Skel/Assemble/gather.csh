#!/bin/csh -f
#foreach f(`sort finished | uniq`)
foreach f(`cat host`)
    echo $f > finished
    rsync -e ssh -avz ${f}:/tmp/kasahara2/"*"-r /tmp/kasahara2/  >& error
end
