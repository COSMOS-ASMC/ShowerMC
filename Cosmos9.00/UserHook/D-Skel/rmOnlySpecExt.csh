#!/bin/csh  -f
#  remove files with only specified file extensions. in the specified dir.
if ( $#argv  <  2  ) then
    echo "Usage: $0 dir {.ext1 .ext2 ...}"  
    echo "where dir is the directory "
    echo ".ext1 .ext2 ... are file name extensions ; the files with those extensions"
    echo "in the directory dir will be removed.
    echo "If no extension is given, nothing happens" 
    exit
endif

set dir = $1
shift

while ( $#argv > 0 )
     rm -f $dir/*$1
     shift
end
