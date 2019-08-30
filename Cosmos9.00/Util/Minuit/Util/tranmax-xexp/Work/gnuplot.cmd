a=0.6246E+13
b=10.49
c=10.98
d=0.8558
x0=562.7
g(x)=a*(x/x0)**b*exp( -c* (x/x0)**d)
xmax = x0* (b/c/d)**(1./d)
#  fix the display region
set xr[xmax-300:xmax+500]
set style arrow 1 size graph 0.02,20 filled linewidth 1.2
set arrow 1 from xmax,g(xmax)/1.9 to xmax,g(xmax) arrowstyle 1
set label 1 sprintf("%6.1f", xmax) center at first xmax,g(xmax)/2
plot "Work/xytemp" ps 2 pt 7,   g(x)  lw 2 lt 3




pause mouse




