set xr[0:]
set log y
set label 1 sprintf("Eg=%10.3e GeV", Eg) at  graph 0.3, graph 0.1 textcolor lt 3
plot "WWW/pair.hist" u 1:2  w his lw 3
pause mouse

#set term png nocrop medium size 480,360
#set term png nocrop medium size 640,480
# next:  can fit 6 graphs in A4 portrate
set term png nocrop medium size  288,216



