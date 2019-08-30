BEGIN {n=0}
$2==name {n++; if(n > howmany) exit;  print $1, $2, $3;next}
n > 0 && n < howmany {print "error: host not enough" > /dev/stderr; exit}
