# to shift (x,y) of +primary by eps (cm)
BEGIN{eps=1.e-5;OFMT="%.12g"; CONVFMT="%.10g"}
NR<3 {print; next}
 {x=$5; y=$6; if(x<0) $5=x+eps; else if(x>0)  $5=x-eps;\
     else {u=rand(); if(u<0.5) $5=eps; else $5=-eps} \
     if(y<0) $6=y+eps; else if(y>0) $6=y-eps; \
     else {u=rand(); if(u<0.5) $6=eps; else $6=-eps} \
     print;
 }
