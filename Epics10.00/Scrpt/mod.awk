BEGIN{nc=0}
NR <= 4 {print;next}
NR>4 && nc==0 {if( $3 == old ) $3=new; nc++; printf("%10d%10d%10d%10d%10d%15.5f\n", $1,$2,$3,$4,$5,$6); next}
NR>4 {print; ++nc;if(nc == 10) nc=0}


