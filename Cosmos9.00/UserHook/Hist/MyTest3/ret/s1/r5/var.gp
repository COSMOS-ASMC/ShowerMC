# variable part of gnuplot commands
plot "E1.dat" using 1:($1**pw*$2) title " s= 0.80 r: 7.943E+00 1.259E+00" w his lw 3,\
  "E3.dat" using 1:($1**pw*$2) title " s= 0.80 r: 7.943E+00 3.162E+00" w his lw 3,\
  "E5.dat" using 1:($1**pw*$2) title " s= 0.80 r: 7.943E+00 7.943E+00" w his lw 3,\
  "E7.dat" using 1:($1**pw*$2) title " s= 0.80 r: 7.943E+00 1.995E+01" w his lw 3
