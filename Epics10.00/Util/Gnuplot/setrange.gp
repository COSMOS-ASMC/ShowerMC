# Mac: getting Work dir is diffcult.
#    Work="`echo /tmp/$USER`"
#    print Work; next NG
#   xmin = `awk '{print $1}' Work/range`

#  use   /tmp/$USER/range directly
xmin = `awk '{print $1}' /tmp/$USER/range`
xmax = `awk '{print $2}' /tmp/$USER/range`
ymin = `awk '{print $3}' /tmp/$USER/range`
ymax = `awk '{print $4}' /tmp/$USER/range`
zmin = `awk '{print $5}' /tmp/$USER/range`
zmax = `awk '{print $6}' /tmp/$USER/range`
set xr[xmin:xmax]
set yr[ymin:ymax]
set zr[zmin:zmax]
#  no gap at the bottom
set ticslevel  0
set term x11 size 600, 750
#  rather horizontal view
set view 85,30,2,1


 
