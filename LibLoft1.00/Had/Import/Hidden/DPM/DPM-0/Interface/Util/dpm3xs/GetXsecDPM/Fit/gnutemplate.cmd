g1(x)=(c1*log(x) + b1) * log(x) + a1
#g2(x)=(c2*log(x) + b2) * log(x) + a2
#g3(x)=(c3*log(x) + b3) * log(x) + a3
#  fix the display region
set log xy
set xr[50:1e11]
plot "Work/xydat" w l lw 2
rep x>50 && x<1.e5 ? g1(x) : 1/0 lw 2
#rep x>1.e5 && x<1.e8 ? g2(x) : 1/0 lw 2
#rep x>1.e8 && x<1.e11 ? g3(x) : 1/0 lw 2

pause mouse

