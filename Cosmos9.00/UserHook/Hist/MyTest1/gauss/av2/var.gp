# variable part of gnuplot commands
plot "1.dat" using 1:($1**pw*$2) title " 4.000E+00 4.000E-01" w his lw 3,\
   "2.dat" using 1:($1**pw*$2) title " 4.000E+00 6.000E-01" w his lw 3,\
   "3.dat" using 1:($1**pw*$2) title " 4.000E+00 8.000E-01" w his lw 3,\
   "4.dat" using 1:($1**pw*$2) title " 4.000E+00 1.000E+00" w his lw 3,\
   "5.dat" using 1:($1**pw*$2) title " 4.000E+00 1.200E+00" w his lw 3
