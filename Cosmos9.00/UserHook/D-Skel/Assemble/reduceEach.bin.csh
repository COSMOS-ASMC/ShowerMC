#!/bin/csh -f
#  Use this script at directory where you have each binary .dat files
#
#  binary  each  .dat is reduced in its size and concatinated.
#  header info is correctly managed: not used now
#  you will get .nrfai-r files for each .dat.
#  you have to combine them by using assemNrfai.csh in Assemble to get
#  reduced .nrfai
#------------fix next
#  .nrfai for all events
setenv NRFAIFILE p18cos0.7.nrfai
#  .dat files
#  setenv HEADER yes
set datfiles=p18cos0.7*.dat
set outfile=/Loft2/Work/kasahara/CosmosData/ForTA/p18New/p18cos0.7/p18cos0.7-all
#-----------------
foreach f($datfiles)
     echo $f >> error
     setenv DATFILE  $f
     (  ~/Cosmos/UserHook/DisParaForTA/Assemble/reduceEachSize.bin.PCLinuxIFC >> $outfile ) >>& error
#     setenv HEADER no
end
