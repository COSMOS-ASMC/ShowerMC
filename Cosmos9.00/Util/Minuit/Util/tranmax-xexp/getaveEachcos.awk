BEGIN{for(i=1; i<=8; i++) sum[i]=0; n=0; cosold=0.}

$9 != cosold && n>0 {print cosold, sum[1]/n, sum[3]/n, sum[5]/n, sum[7]/n,  n;  \
		     for(i=1; i<=8; i++) sum[i]=0; n=0};
{for(i=1; i<=8; i++) sum[i]+=$i; cosold=$9; n++}
END  {if(n>0) print cosold, sum[1]/n, sum[3]/n, sum[5]/n, sum[7]/n,  n}



