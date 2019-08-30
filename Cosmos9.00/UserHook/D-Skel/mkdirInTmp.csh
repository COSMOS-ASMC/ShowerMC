#!/bin/csh -f
#  make /tmp/$USER for every host
if ( $#argv == 0 ) then
    echo "./mkdirInTmp.csh hostfile"
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
    ssh $f mkdir -p /tmp/$USER
end
echo "Number of cpu = $n"
    
