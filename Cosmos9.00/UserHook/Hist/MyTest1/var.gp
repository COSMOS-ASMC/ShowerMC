# variable part of gnuplot commands
plot "1.dat" using 1:($1**pw*$2) title "" w his lw 3,\
   "2.dat" using 1:($1**pw*$2) title "" w his lw 3,\
   "3.dat" using 1:($1**pw*$2) title "" w his lw 3
