BEGIN { s = 0.5; t = 0. }
NR > 18 && NF == 10 { 
    if(s > 1.5) s= 0.5;
    print $1, $2, $3, $4, $5, $6, $7,  s, t,0; 
    s = s + 1.0/36.0;next
  }
NF==0 && NR>18 { print; t = t+ 1./18.0;next }
{print}
