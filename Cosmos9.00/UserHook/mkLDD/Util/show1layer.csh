#!/bin/csh
#  to show 1 line of .hyb data  at specified depth
if ( $#argv != 1 ) then
    echo "Usage: ./show1layer.csh depth"
    exit
endif
foreach f(*.hyb)
   echo $f
   awk -f $COSMOSTOP/UserHook/mkLDD/Util/show1layer.awk depth=$1  $f
end
