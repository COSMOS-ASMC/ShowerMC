NF==0 && n==2 {print tetacm, tetalab; n=0; next};
NF==0 {n=0;next};
$1==1 && $2==6 && $4==0 \
{tetacm=$10; tetalab=$11; next} 
#  {pt=sqrt($5**2+$6**2); teta=atan2(pt,$7)*180./3.14;n++;next}; 
{n++}
