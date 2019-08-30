#  this is almost the same as geomvewconfi2.awk
#  but  only for non-subdetectors
BEGIN {nonsubd  = -1; subd = -1}
NR == 1{next}
NF > 2 && ( $1 == "CMESH" || $1 == "OFF") && $6==0  {
  nonsubd = 1;
  subd = 0;
  nc = $4;
  media = $5;
  if(nc < 10 ) seq = "0000"nc;
  if(nc >=10 && nc < 100) seq ="000"nc;
  if(nc>=100 && nc < 1000) seq ="00"nc;
  if(nc>=1000 && nc < 10000) seq ="0"nc;
  if(nc>=10000) seq = nc;
  print "LIST" > outdir"/"media"_"seq".list";
}

NF > 2 && ( $1 == "CMESH" || $1 == "OFF") && $6>0  {
  subd = 1;
  nonsubd = 0;
  file = $7$6".list"

  if( system( "echo test -f " outdir"/"file " | sh " ) != 0 )  print "LIST" > outdir"/"file
}

nonsubd == 1  {print > outdir"/"media"_"seq".list"; next}
subd == 1 {print >> outdir"/"file}
     

