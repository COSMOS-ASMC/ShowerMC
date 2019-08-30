#!/bin/csh
if( $#argv != 1)  then
   echo "Usage: ppmtojpeg dir_for_ppmfiles"
   exit
endif
cd $1;
    foreach  f(ts*.ppm)
           set g = $f:r
           ppmtojpeg $f > $g.jpeg
    end

			       
