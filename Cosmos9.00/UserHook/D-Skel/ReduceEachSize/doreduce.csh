#!/bin/csh -f
setenv NRFAIFILE ../Assemble/p19cos0.9.nrfai
foreach f(../Assemble/OutDir/*.dat)
   setenv DATFILE $f
  ./reducebinSizePCLinuxIFC  > $f.reduced
end
