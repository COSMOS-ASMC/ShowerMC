#!/bin/csh

set  head=${COSMOSTOP}/UserHook/mkLDD/Util

if ( $#argv != 1 ) then
cat <<EOF 
   Usage:  ~/Cosmos/UserHook/mkLDD/Util/showMinMaxAge.csh depth
   To get min and max age and cog values  at a give depth
   among the .hyb files in the current direcotry
   and print them with the file name
EOF
  exit 1
endif
if( -f temp$$ ) then
    rm -f temp$$
endif



foreach f(*.hyb)
   awk -f $head/show1layer.awk depth=$1  $f | awk '{print $5, $6, file}'   file=$f  >>  temp$$
end

awk -f $head/getMinMaxAge.awk i=1  temp$$
awk -f $head/getMinMaxAge.awk i=2  temp$$

rm -f temp$$
