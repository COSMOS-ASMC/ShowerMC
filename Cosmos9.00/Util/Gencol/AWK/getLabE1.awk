$1~/\$HPARAM/, $1~/\$END/ {if($1~/UserHookc/) {
  l1=index($0,"'");  d=substr($0,l1+1); l2=index(d,"'");
  dd1=substr(d,1,l2-1); print dd1;
  found=1; next;}
 else if(found==1 ) {
  l1=index($0,"'");  d=substr($0,l1+1); l2=index(d,"'");
  dd2=substr(d,1,l2-1); print dd2;
 }
}


