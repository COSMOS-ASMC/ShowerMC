set xr[0:1]
unset key
set label 1 sprintf("%s:From Eg=%4.3e GeV; log10 step=%4.3f",media, Eg,step) at  graph 0.1, graph 0.95 textcolor lt 3
set xlab "k=Ek/(Eg-2me)" off 19,0.3
set ylab "ds/dk"   off 0.3,6 
plot  "WWW/fpair.func" w l lw 3 lt 3
pause mouse


set term png nocrop medium size 480,360
#set term png nocrop medium size 640,480
# next:  can fit 6 graphs in A4 portrate
# set term png nocrop medium size  288,216






