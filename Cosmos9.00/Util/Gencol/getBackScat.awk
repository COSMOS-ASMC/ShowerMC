$1=="h" && NR==1 {print;  next};
$1=="h" {print "s ", sumG, bkpiK, bkp;
    sumG=0; bkpiK=0; bkp=0; print; next}
{bdeg=atan2($4,$8)*180/3.141592;
    if( $1==4 && $3==0) sumG += $7;
    else if( $1==1 ) sumG +=  $7;
    if(bdeg >90 && $3 !=0 ) {
	if($1==4 || $1==5 )  bkpiK++;
	else if( $1==6 )  bkp++;
	print $1, $2,$3, $7, bdeg, $8;
    }
}
