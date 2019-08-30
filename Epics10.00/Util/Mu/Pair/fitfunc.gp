#!/usr/local/bin/gnuplot -persist
#
#    	plot window:      hit 'h'
xmin(E)=4*0.511/E
x0(E)=25*0.511/E
x1(E)=0.01/(E/3000)**0.4
x2(E) = 0.05/(E/3000)**0.16
x3(E) =0.8
f0(x) = a0*( x/xmin(E) -1.)**b0 * (x/x0(E))**(-c0)
f1(x)=a1*(x/x0(E))**(-b1)* (1+x/sqrt(x1(E)*x0(E)))**(-c1)
#f1(x)=a1*(x/x0(E))**(-b1)* (1+(x/sqrt(x2(E)*x1(E))))**(-c1)
#f2(x)=a2*(x/x1(E))**(-b2)* (1+(x/x2(E)))**(-c2)
f2(x)=a2*(x/x1(E))**(-b2)* (1+(x/sqrt(x1(E)*x2(E))))**(-c2)
f3(x)=a3*(x/x2(E))**(-b3)* (1+x/sqrt(x2(E)*x3(E)))**(-c3)
#    EOF
