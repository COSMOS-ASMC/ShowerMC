# variable part of gnuplot commands
plot "1.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 1.259E+00" w his lw 3,\
plot "E1.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 1.259E+00" w his lw 3,\
   "2.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 3.162E+00" w his lw 3,\
  "E3.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 3.162E+00" w his lw 3,\
   "3.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 7.943E+00" w his lw 3,\
  "E5.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 7.943E+00" w his lw 3,\
   "4.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 1.995E+01" w his lw 3,\
  "E7.dat" using 1:($1**pw*$2) title " s= 1.20 r: 5.012E+01 1.995E+01" w his lw 3
