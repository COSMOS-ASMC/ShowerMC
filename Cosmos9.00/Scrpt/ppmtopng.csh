#!/bin/csh
if( $#argv == 0 ) then
    echo "Usage: ppmtopng dir_for_skelfile [compression(0-9)]"
     exit
endif
    cd $1;
    foreach  f(ts*.ppm)
       set g = $f:r
      if( $#argv == 2 ) then
         pnmtopng -comp $2  $f > $g.png
      else
         pnmtopng   $f > $g.png
      endif
    end
