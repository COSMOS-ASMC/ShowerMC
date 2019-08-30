#  this is to calculate rho(r)*2pir at 6(=nr) r's
#  by using .lat data. We assume age and cog are given
#  at the first line. 
#  The ouput format will be
#    code fai age cog  r,y, r,y... r,y  (6 pairs; y= 2pixrxrho(r))
#     


BEGIN{ra[1]=5.0; ra[2]=7.0; ra[3]=10.0; ra[4]=20.0; ra[5]=30.0; ra[6]=50.0; nr=6 }
NR==1 {age=$1; cog=$2; next}
#$1=="l" {code=$3; fai=$4; lc=0; if(code == 3) pw=pwin; else pw=0.5; next }
$1=="l" {code=$3; fai=$4; lc=0;next }
lc==0 {lc++;next}
lc==1 {lc++;next}
lc==2 {lc++;a=$1; b=$2; c=$3; pw=$4;
       for(i=1; i<=2; i++) { r=ra[i]; rho[i]=a/(r**(b+c*r**pw))}; next}
lc==3 {lc++; a=$1; b=$2; c=$3; pw=$4;
       for(i=3; i<=nr; i++) { r=ra[i]; rho[i]=a/(r**(b+c*r**pw))}
       printf("%d %d %6.3f %6.3f", code, fai, age, cog);
       for(i=1; i<=nr; i++) {printf(" %6.1f  %11.3e", ra[i], rho[i])}
       print " "
}
