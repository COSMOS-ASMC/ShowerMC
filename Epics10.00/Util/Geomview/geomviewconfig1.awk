# input   data format
#    x y z      cn  struc    subdidx subdname
#    x y z      cn  struc    subdidx subdname
#
#    x y z      cn  struc    subdidx subdname
#    x y z      cn  struc    subdidx subdname
#
#
#    x y z      cn  struc    subdidx subdname
#    x y z      cn  struc    subdidx subdname
#
#    x y z      cn  struc    subdidx subdname
#    x y z      cn  struc    subdidx subdname
#    
#  strategy.  store  the data in an array for the same compnent.
#             if struc is box, treat specially,
#             otherwsie convrt them into CMESH
# output; CMESH data and/or COFF data are on stdout
#
#   if complete box exists, COOF data is put as
#
#    OFF #   struc compn media subdidx subdname
#    8  6  12
#    lines for the box
#
#    CMESH #  struc compn media subdidx subdname
#    lines for CMESH
# ...
#
BEGIN {first= 0; lc = 0; struc = " "; compn = " " ; subdidx = 0; subdname=" " }
NF==0 {lc++; lines[lc] = $0; next}

$4 != compn && first > 0 {
  save = $0;
#      if full box, treat it specially 
  if( struc == "box" && lc > 42) procbox(lines); 
  else if( substr(struc,1,8) == "fpolygon" ){
      if( concave == "ok")  procfpolygon(lines,lc);     procother(lines, lc);
  }
  else {
    procother(lines, lc);
  }
  lc = 0;
  $0 =save;
}

first == 0 {first = 1 }

{lc++;  lines[lc] = $0 ; compn=$4; struc = $5; subdidx = $6; subdname=$7;next}


END {if(lc > 42 && struc == "box") procbox(lines);
    else if( substr(struc,1,8) == "fpolygon" ) {
	if( concave == "ok") procfpolygon(lines,lc);       procother(lines, lc);
    }
    else {
       procother(lines, lc);
     }
   }

#  ****************************************
function procother(lines, lc) {
#   all mesh  other  than complete box

if(lc <= l) return;
blc = 0; prev="0"; u=0; v=1; yet=0;
for(i=1; i<=lc; i++) {

  $0=lines[i];


  if( NF == 0 ) { blc++; prev="0"; continue}
  if(prev=="0" && blc==1) { blc=0;  v++; u=0; prev="d"}

  if(prev=="0" && blc >= 2 && u>0 ){
    blc=0;
    if(yet == 0 ) {
      print "CMESH  #  ", struc, compn, media, subdidx, subdname;
      yet = 1
    }
    else
      print "CMESH";

    print u, v;
    for(j=1; j<=v; j++){
      for(k=1; k<=u; k++) print x[k,j], y[k,j], z[k,j], c1, c2, c3, c4;
      print " ";
    };
    u=0; v=1; prev="d"
  }
  u++; x[u,v]=$1; y[u,v]=$2; z[u,v]=$3;


#  u++; x[u,v]=$1; y[u,v]=$2; z[u,v]=$3; compn=$4; struc=$5; subdidx=$6; subdname=$7;
}
if(u>0) {


  if(yet == 0) 
    print "CMESH # ", struc, compn, media,  subdidx, subdname;
  else
    print "CMESH"

  print u, v;
  for(j=1; j<=v; j++){
    for(i=1; i<=u; i++) print x[i,j], y[i,j], z[i,j], c1, c2, c3, c4;
    print " "
  };
}
}
#  *******************************
function procbox(lines) {
#  complete box
#  print "LIST"
  print "OFF #  ", struc, compn, media, subdidx, subdname;
  print "8 6 12";
#      output 8 vertex
  $0 = lines[1];
  print $1, $2, $3;
  $0 = lines[2];
  print $1, $2, $3;

  $0 = lines[4];
  print $1, $2, $3;
  $0 = lines[5];
  print $1, $2, $3;

  $0 = lines[11];
  print $1, $2, $3;
  $0 = lines[12];
  print $1, $2, $3;
  
  $0 = lines[39];
  print $1, $2, $3;
  $0 = lines[40];
  print $1, $2, $3;
#  output face info. with color index
#  verteces    r g b
  print "4 3 1 0 2 " c1, c2, c3, c4
  print "4 4 5 7 6 " c1, c2, c3, c4
  print "4 2 3 7 6 " c1, c2, c3, c4
  print "4 0 1 5 4 " c1, c2, c3, c4
  print "4 0 4 6 2 " c1, c2, c3, c4
  print "4 1 5 7 3 " c1, c2, c3, c4
}
function procfpolygon(lines,lc) {
#  flat polygon
  print "OFF #  ", struc, compn, media, subdidx, subdname;
  nv=(lc-4);
  print nv,  2,   nv;
#     
  for(i= 1; i<= lc;  i++) {
      $0 = lines[i];
      if( $0 != " " )  print $1, $2, $3;
  }
#  output face info. with color index
#  verteces    r g b
  printf("%d ", nv/2);
  for(i=1; i<=nv/2; i++) {
      printf("%d ", i-1);
  }    
  printf("  %d %d %d %d\n", c1, c2, c3, c4);

  printf("%d ", nv/2);
  for(i=nv/2+1; i<=nv; i++) {
      printf("%d ", i-1);
  }    
  printf("  %d %d %d %d\n", c1, c2, c3, c4);
}
