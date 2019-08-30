#!/bin/csh
foreach f(*.dat)
(  echo $f |  ../bin2asciiPCLinuxIFC > ../Ascii/$f ) >> & error
end
