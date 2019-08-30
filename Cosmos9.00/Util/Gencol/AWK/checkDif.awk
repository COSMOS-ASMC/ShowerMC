$1=="h" {if($3==1) {sd=1;next} \
         if( $4==1) {sd=2;next};    \
         sd=0; next}
$1=="p" && sd!=0 {print sd,  $0; if( $3!=0) print "****"}



