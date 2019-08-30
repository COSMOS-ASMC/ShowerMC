# set log y
set xr[1.e-5:1]
unset key
set label 1 sprintf("%s:From Ek=%4.3e GeV; log10 step=%4.3f",media, Ek,step) at  graph 0.1, graph 0.95 textcolor lt 3
set xlab "k=Eg/Ek"
set ylab "kds/dk"
plot  "Work/fbrems.func" w l lw 3 lt 3
pause mouse


set term png nocrop medium size 480,360
#set term png nocrop medium size 640,480
# next:  can fit 6 graphs in A4 portrate
# set term png nocrop medium size  288,216






