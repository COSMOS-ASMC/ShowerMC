BEGIN{sum=0; n=0; cosold=0.}

$3 != cosold && n>0 {print sum/n, cosold, n; sum=0; n=0};
 {sum+=$1; cosold=$3; n++}
END  {if(n>0) print sum/n, cosold, n}


