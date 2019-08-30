#!/bin/csh
    foreach  f(ts*.ppm)
       echo $f
       set g = $f:r
       pnmquant 256 $f |  ppmtogif  > $g.gif
    end
