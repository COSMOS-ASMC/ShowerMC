set xr[0:1]
set yr[0:]
unset key
set label 1 sprintf("%s:From Eg=%4.3e GeV; log10 step=%4.3f",media, Eg,step) at  graph 0.1, graph 0.95 textcolor lt 3
set xlab "k=Ek/(Eg-2me)"  off 19,0.3
set ylab "ds/dk"   off   0.5,6
