#!/bin/csh  -f
@ n=1
foreach f(*.hyb)
     echo $f
     set nn=`echo $n | awk '{printf("%3.3d\n",$1)}'`
     echo $nn
     echo $f:r:r.$nn.hyb
     mv $f  $f:r:r.$nn.hyb
     @ n++
end
