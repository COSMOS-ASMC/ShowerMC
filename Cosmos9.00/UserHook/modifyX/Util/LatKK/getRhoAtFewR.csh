#!/bin/tcsh
if ( -f rhoAtRs.data ) then
  rm -f rhoAtRs.data
endif
# set pwa=(`grep "%pw%" ~/Cosmos/UserHook/Minuit/Util/latFit/latFit.f`)
# set pw=$pwa[3]
# echo "pw for muon in R> 10 is $pw"

foreach f($1/*.lat)
#     awk -f getRhoAtFewR.awk pwin=$pw $f >>  rhoAtRs.data
     awk -f getRhoAtFewR.awk  $f >>  rhoAtRs.data
end







