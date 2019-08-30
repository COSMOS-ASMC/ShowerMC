#!/bin/csh
if( $#argv == 0 )  then
    echo "Usage: ./ppmtogif.csh ppmfile_dir [s]"
    echo "optional s is to skip color map making for speed-up "
    echo "Use s, if you are sure that all ppm files contain < 256 colors"
    exit
endif
cd $1;
    foreach  f(ts*.ppm)
       echo $f
       set g = $f:r
       if( $2 == "s" ) then
            ppmtogif $f  > $g.gif
       else
            pnmquant 256 $f |  ppmtogif  > $g.gif
       endif	    
    end
