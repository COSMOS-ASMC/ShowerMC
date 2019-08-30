BEGIN {nevc=0}
$1==code && $2==fai && $3>minAge && $3<maxAge  \
{for(i=5;i<=16;i+=2) print $i, $(i+1)} 

