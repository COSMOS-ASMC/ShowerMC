NR == 1 {next}
NF > 2 && ( $1 == "CMESH" || $1 == "OFF")  {
  nc = $4;
  media = $5;
  if(nc < 10 ) seq = "0000"nc;
  if(nc >=10 && nc < 100) seq ="000"nc;
  if(nc>=100 && nc < 1000) seq ="00"nc;
  if(nc>=1000 && nc < 10000) seq ="0"nc;
  if(nc>=10000) seq = nc;
  print "LIST" > outdir"/"media"_"seq".list";
}

{print > outdir"/"media"_"seq".list"}





