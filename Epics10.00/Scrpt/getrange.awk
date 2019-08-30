BEGIN { xmin=1.e5; xmax=-1.e5; ymin=1.e5; ymax=-1.e5; zmin=1.e5;
    zmax=-1.e5;}
NF>0 {
  if($1 < xmin) xmin=$1;	
  if($1> xmax)  xmax =$1; 
  if($2 < ymin) ymin=$2; 
  if($2> ymax)  ymax =$2;
  if($3 < zmin)  zmin = $3;
  if($3 > zmax)  zmax = $3;
  Dx = (xmax - xmin); Dy=(ymax - ymin);
  if( Dx > Dy ) { 
     h= (Dx-Dy)/2; ymax=ymax+h; ymin =ymin-h; 
  } 
  else { 
     h=(Dy-Dx)/2; xmax=xmax+h; xmin =xmin-h;
#     print "h, xmax, xmin",h, xmax, xmin; 
  }; 
  print xmin, xmax, ymin, ymax, zmin, zmax > file;
  exit;
}

