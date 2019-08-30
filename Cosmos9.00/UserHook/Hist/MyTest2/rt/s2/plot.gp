# gnuplot commands
set title " Time spectrum" font "Times,20"
set xlabel "T(ns)" offset  0,0.5 font "Times-Italic,22"
pw ="0.00"
set ylabel "dN/NdT(/ns)" offset  2.5,0  font "Times-Italic,22"
set grid mxtic xtic  mytic ytic
set log x
set log y
#  set key x,y
#  set xrange [1:10]
#  set yrange [1:10]
# tab is: y dn/dx/dy dn/dx/dy* dn  dy x
# tab is: z dn/dx/dy/z dn/dx/dy/dz* dn dz x y
# graph(dn/dx..) is normalzied.
# To show unnormalzied one(dn/dx..*); change $2-->$3 in the last call command
# To change the line style etc, modify the last "w his"
call "var.gp"
