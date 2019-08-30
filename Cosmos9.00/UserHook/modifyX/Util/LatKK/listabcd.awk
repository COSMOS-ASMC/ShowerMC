#BEGIN{q="'"; print q"lateral fit"q; 
# print q"code"q q"region"q q"a"q q"b"q q"c"q q"pw"q q"age"q q"cog"q q"maxdiff"q}

NF==2 {age=$1; cog=$2; ok=0; next}
$1=="l" &&  $4==fai {rgc=0; code=$3;ok=1;next}
ok==1  {rgc++; a=$1; b=$2; c=$3; pw=$4; maxdiff=$5;
	den=0.;
	if(R < 0.1 && rgc==1) den=a/( R**(b+c*R**pw));
	if(R < 0.1  && R>0.01 &&  rgc==2) den=a/( R**(b+c*R**pw));
	if(R < 1 &&  R>0.1 && rgc==3) den=a/( R**(b+c*R**pw));
	if(R> 1. && rgc==4)  den=a/( R**(b+c*R**pw));
	print code, rgc,$1,$2,$3,$4, age, cog, $5, R, den;
	if(rgc==4) ok=0; next
}



