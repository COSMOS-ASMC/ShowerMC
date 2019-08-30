NF>0 {
  xmin=-5e3; if($1 < xmin) xmin=$1;	
  xmax=5e3;  if($1> xmax)  xmax =$1; 
  ymin=-5e3; if($2 < ymin) ymin=$2; 
  ymax=5e3;  if($2> ymax)  ymax =$2;
  zmax =$3;  zmin=0; 
  Dx = (xmax - xmin); Dy=(ymax - ymin);
  if( Dx > Dy ) { 
     h= (Dx-Dy)/2; ymax=ymax+h; ymin =ymin-h; 
  } 
  else { 
     h=(Dy-Dx)/2; xmax=xmax+h; xmin =xmin-h;
     print "h, xmax, xmin",h, xmax, xmin; 
  }; 
  print xmin, xmax, ymin, ymax, zmin, zmax > file;
  exit;
}

