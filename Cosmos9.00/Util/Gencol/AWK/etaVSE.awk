# tx = sint sinf 
# ty = sint cosf
# tx = cost 
 
NR>2 && NF> 0 {E=$4; code=$1; subc=$2; tx=$5; ty=$6; cost=$7; sint = sqrt(tx**2 + ty**2);
    teta = atan2(sint,cost);  if(code==9) Epn=E/subc; if(code != 9)  Epn=E;
    hteta= teta/2;      eta = - log( sin(hteta)/cos(hteta) );
    print  eta, Epn
}

