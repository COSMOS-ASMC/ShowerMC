#!/bin/csh
if( $#argv == 0 ) then
    echo "Usage: ppmtopng dir_for_skelfile"
     exit
endif
    cd $1;
    foreach  f(ts*.ppm)
       set g = $f:r
       pamtotga   $f > $g.tga
    end
