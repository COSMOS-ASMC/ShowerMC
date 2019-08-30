$1=="h" {nev++; layer=0}
$1=="t" { if(nev>1) { if($3 != gram[$2]) {print -1000 > /dev/stderr;  exit};}
	  gram[$2]=$3; \
	  sumg[$2]+=$7; sume[$2]+=$8; \
	  summu[$2]+=$9; sumh[$2]+=$10; \
	  sumhyb[$2]+=$11; sumdEdx[$2]+=$12; sumdE[$2]+=$13; \
	  sumDead[$2]+=$14; if(layer < $2) layer=$2}
END {for(i=1;i<=layer; i++) print i, gram[i], sumg[i]/nev, sume[i]/nev, 
     summu[i]/nev, sumh[i]/nev, sumhyb[i]/nev, 
     sumdEdx[i]/nev,  sumdE[i]/nev, sumDead[i]/nev
}
