{i++; x[i]=$1; y[i]=$2};
END {for(j=1;j<=i;j++) { printf("%5.3f, ", x[j]); 
	if( j == int(j/5)*5 ) print " ";}
    print " ";
    for(j=1;j<=i;j++) { printf("%5.3f, ", y[j]);
	if( j == int(j/5)*5 ) print " ";} print " "
}
