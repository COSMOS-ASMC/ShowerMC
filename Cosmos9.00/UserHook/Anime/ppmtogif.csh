#!/bin/csh
    foreach  f(ts*.ppm)
       set g = $f:r
       ppmtogif $f > $g.gif
    end
