# variable part of gnuplot commands
plot "1.dat" using 1:($1**pw*$2) title " s= 1.00 r: 7.943E+00 1.259E+00" w his lw 3,\
   "2.dat" using 1:($1**pw*$2) title " s= 1.00 r: 7.943E+00 3.162E+00" w his lw 3,\
   "3.dat" using 1:($1**pw*$2) title " s= 1.00 r: 7.943E+00 7.943E+00" w his lw 3,\
   "4.dat" using 1:($1**pw*$2) title " s= 1.00 r: 7.943E+00 1.995E+01" w his lw 3
