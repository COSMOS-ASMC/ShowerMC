#!/bin/csh
#  to get min and max age and cog values among the  .hyb files in the
#  current direcotry and print them with the file name
#  Usage:  $COSMOSTOP/UserHook/mkLDD/showMinMaxAge.csh depth
# at a give depth

if ( $#argv != 1 ) then
    echo "Usage: ./showMinMaxAge.csh depth"
    exit
endif
if( -f temp$$ ) then
    rm -f temp$$
endif

foreach f(*.hyb)
   awk -f $COSMOSTOP/UserHook/mkLDD/show1layer.awk depth=$1  $f | awk '{print $5, $6, file}'   file=$f  >>  temp$$
end

awk -f $COSMOSTOP/UserHook/mkLDD/getMinMaxAge.awk i=1  temp$$
awk -f $COSMOSTOP/UserHook/mkLDD/getMinMaxAge.awk i=2  temp$$

rm -f temp$$




