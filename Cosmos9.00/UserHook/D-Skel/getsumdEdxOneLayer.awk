# to  get sum of  number of ptcls of a given type  at a gvien layer from .nrfai 
#   awk -f getsumdEOneLayer.awk   layer=yy xxx.nrfai
#
BEGIN {sum=0; ok=-1}
$1=="dE/dx" && $3==layer  { ok=1; next}
ok==1 {for(i=1; i<=NF; i++) sum+=$i; if(NF==2) ok=0; next}
sum > 0 {print sum; exit}


