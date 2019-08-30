# variable part of gnuplot commands
plot "r1.dat" using 1:($1**pw*$2) title " s= 0.80 r: 1.259E+00" w his lw 3,\
  "r3.dat" using 1:($1**pw*$2) title " s= 0.80 r: 3.162E+00" w his lw 3,\
  "r5.dat" using 1:($1**pw*$2) title " s= 0.80 r: 7.943E+00" w his lw 3,\
  "r7.dat" using 1:($1**pw*$2) title " s= 0.80 r: 1.995E+01" w his lw 3,\
  "r9.dat" using 1:($1**pw*$2) title " s= 0.80 r: 5.012E+01" w his lw 3
