$1=="h" {next}
$8 > 10 {teta=atan2($4,$8)*180/3.1415; print $7, teta, $1,$3}
