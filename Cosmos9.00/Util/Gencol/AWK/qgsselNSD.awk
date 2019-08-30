# for qgsjet1 and 2
BEGIN {s=0} 
$1=="h" && $2==1  {s=1;nev++; next}
NF>0 &&  s==1  && $2!=0 {print $4;next}
NF==0 { s=0 }
END {print nev >> "/tmp/kasahara/nev"}
