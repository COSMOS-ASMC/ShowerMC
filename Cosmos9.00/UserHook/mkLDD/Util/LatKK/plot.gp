norm1=322
norm2=401
second=2

set term aqua fs 24

set log xy
plot "gammaS0-1.data" ps 0.6
replot "../../T1-7/V4/gammaS0-1.data" u ($1*1.03):2 ps 0.6 lt 3

pause second

plot "elecS0-1.data" ps 0.6
replot "../../T1-7/V4/elecS0-1.data" u ($1*1.03):2 ps 0.6 lt 3

pause second

plot "muonS0-1.data" ps 0.6
replot "../../T1-7/V4/muonS0-1.data" u ($1*1.03):2 ps 0.6 lt 3
pause second

plot "hadronS0-1.data" ps 0.6
replot "../../T1-7/V4/hadronS0-1.data" u ($1*1.03):2 ps 0.6 lt 3

pause second

unset log xy
plot "gammaS0-1at50.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/gammaS0-1at50.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second

plot "gammaS0-1at30.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/gammaS0-1at30.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second

plot "gammaS0-1at20.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/gammaS0-1at20.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second

plot "gammaS0-1at10.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/gammaS0-1at10.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second

plot "elecS0-1at50.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/elecS0-1at50.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second

plot "elecS0-1at30.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/elecS0-1at30.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second

plot "elecS0-1at20.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/elecS0-1at20.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second

plot "elecS0-1at10.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/elecS0-1at10.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second


plot "muonS0-1at50.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/muonS0-1at50.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second

plot "muonS0-1at30.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/muonS0-1at30.hist" u 1:($2/norm1)  w his lw 3 lt 3

pause second


plot "muonS0-1at20.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/muonS0-1at20.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second

 plot "muonS0-1at10.hist" u 1:($2/norm2) w  his lw 3
 replot "../../T1-7/V4/muonS0-1at10.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second


plot "hadronS0-1at50.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/hadronS0-1at50.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second

plot "hadronS0-1at30.hist" u 1:($2/norm2) w  his lw 3
 replot "../../T1-7/V4/hadronS0-1at30.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second

plot "hadronS0-1at20.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/hadronS0-1at20.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second

plot "hadronS0-1at10.hist" u 1:($2/norm2) w  his lw 3
replot "../../T1-7/V4/hadronS0-1at10.hist" u 1:($2/norm1)  w his lw 3 lt 3
pause second

