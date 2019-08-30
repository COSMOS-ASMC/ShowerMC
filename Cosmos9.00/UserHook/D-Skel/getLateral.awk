#   extract some dato of a fixed layer data and use 
#  awk -f getLateral.awk data > newdata

BEGIN {r=0.01};
{for(i=1;i<=NF;i++) { print r,  $i; r=r*10.**0.1}; if(NF==2) r=0.01}

