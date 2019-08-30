#!/bin/csh
    foreach  f(ts*.ppm)
       set g = $f:r
       ppmtojpeg $f > $g.jpeg
    end
