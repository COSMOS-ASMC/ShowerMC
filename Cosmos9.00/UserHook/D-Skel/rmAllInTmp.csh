#!/bin/csh -f
if ( $#argv != 2 ) then
    echo "./rmAllInTmp.csh hostfile  dir"
    echo "All files in dir will be deleted"
    echo "If dir is absent, it will be created"
    exit
endif
set fold="xx"
@ n = 0
foreach f(`awk '{print $2}' $1`)
    if ( "x$f" == "x#" ) continue
    if ( "x$f" == "x" ) continue
    @ n++;
    if ( $f == $fold ) continue
    set fold=$f
    echo $f
    ssh $f mkdir -p $2
    ssh $f rm -f $2/"*"
end
echo "Number of cpu = $n"
