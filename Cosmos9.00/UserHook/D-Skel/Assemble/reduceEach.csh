#!/bin/csh -f
#  Use this script at directory where you have each .dat files
#
#  ascii each  .dat is reduced in its size and concatinated
#  header info is correctly managed; not  used now
#  you will get .nrfai-r files for each .dat.
#  you have to combine them by using assemNrfai.csh in Assemble to get
#  reduced .nrfai
#------------fix next
#  .nrfai for all events
setenv NRFAIFILE p20cos0.95.nrfai
#  .dat files
set datfiles=p18cos0.9*.dat
set outfile=p18cos0.9-all
#-----------------
foreach f($datfiles)
#  set exist=`grep  $f datfilelist`
#  if( x$exist == "x" ) then
#     echo $f >> datfilelist
     echo $f >> error
     setenv DATFILE  $f
     (  ~/Cosmos/UserHook/DisParaForTA/Assemble/reduceEachSizePCLinuxIFC >> $outfile ) >>& error
#  endif
end
