# gnuplot commands
set title " Test Gaussian dist." font "Times,20"
set xlabel "x(m)" offset  0,0.5 font "Times-Italic,22"
pw ="0.00"
set ylabel "dN/Ndx(/m)" offset  2.5,0  font "Times-Italic,22"
set grid mxtic xtic  mytic ytic
#  set key x,y
#  set xrange [1:10]
#  set yrange [1:10]
# tab is: x dn/dx dn/dx* dn  n(>x) <x> dx
# graph(dn/dx..) is normalzied.
# To show unnormalzied one(dn/dx..*); change $2-->$3 in the last call command
# To change the line style etc, modify the last "w his"
call "var.gp"
