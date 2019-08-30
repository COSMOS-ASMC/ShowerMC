function outptcl(file, x,y,z,n, code, chg) {
    print "SKEL" > file ;
    print  n,  n > file ;
    for(i=1; i<=n; i++){ print x[i], y[i], z[i] >> file}; 
    for(i=1; i<=n; i++) { print 1, i-1, 1, 0, 0 >> file };

 }

$4==2 && $5 < 0 {ne++; xne[ne]=$1;yne[ne]=$2;zne[ne]=$3;next}
$4==2 && $5 > 0 {pe++; xpe[pe]=$1;ype[pe]=$2;zpe[pe]=$3;next}
$4==3 && $5 < 0 {nm++; xnm[nm]=$1;ynm[nm]=$2;znm[nm]=$3;next}
$4==3 && $5 > 0 {pm++; xpm[pm]=$1;ypm[pm]=$2;zpm[pm]=$3;next}
$4  >3    {o++; xo[o]=$1; yo[o]=$2; zo[o]=$3;next}
END {
  if(ne >0 ) outptcl("e-.dat", xne, yne, zne, ne,2, -1);
  if(pe >0 ) outptcl("e+.dat", xpe, ype, zpe, pe,2, 1);
  if(nm >0 ) outptcl("mu-.dat", xnm, ynm, znm, nm, 3, -1);
  if(pm >0 ) outptcl("mu+.dat", xpm, ypm, zpm, pm, 3, 1);
  if( o >0 ) outptcl("other.dat", xo, yo,  zo, o, 0, 0);
}
