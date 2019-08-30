# set log y
unset key
set xr[1.e-4:1]
set format y "%.2e"
set xlabel "k=Eo/KEmu" off 19,0.3
set label 1 sprintf("muon nuclear int. in %s,\n Ek=%12.3e GeV #=%i", media,Ek,nev) at  graph 0.04, graph 0.1 textcolor lt 3

plot "WWW/nuci.hist" u 1:($2*$1)  w his lw 3, "WWW/nuci.func" w l lw 3 lt 3
pause mouse

#set term png nocrop medium size 480,360
#set term png nocrop medium size 640,480
# next:  can fit 6 graphs in A4 portrate
set term png nocrop medium size  288,216



