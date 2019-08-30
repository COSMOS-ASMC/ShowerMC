#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#define  NMAX  32768
/*     usage:
 *       histo [-l] xmin bin norm
 * standard input is assumed.
 * This creates a table of 
 *
 *   xcenter  dN/dx/norm  dN  N(>x)  <x>
 *
 *   where xcenter is  xmin+bin/2 xmin+bin+bin/2, ... 
 *   dN is the number of data in [xcenter-bin/2, xcenter+bin/2)
 *   N(>x) is the number of data > x-bin/2.
 *   if -l is usued, all counting is done taking log10
 *   of the variable values but bin is understood as
 *   the one in log10 scale (say, if  0.1, it  means
 *   log10(x) is binned by 0.1 step )
 *   if norm is missing, 1 is assumed.
 */
int main(int argc, char *argv[])
{
  int lg;
  float xmin, bin, to10, x, xm, norm, xx;
  long int sum, nmax, imax, imin, i;
  float inc, dx;
  long int dn[NMAX];
  float mean[NMAX];


  if(argc < 3) {
    fprintf(stderr,"histo [-l] min bin [norm] < infile > outfile\n");
    exit(1); 	
  }
  lg =  strcmp(argv[1], "-l") ? 0 : 1 ;
  if(argc < 3+lg || argc > 4+lg ) {
    fprintf(stderr,"histo [-l] min bin [norm] < infile > outfile\n");
    exit(1); 	
  }
  nmax = NMAX - 1;
  xmin = atof(argv[lg+1]);
  bin = atof(argv[lg+2]);


  norm = (argc == 3+lg) ? 1.0 : atof(argv[lg+3]);
/*   
  fprintf(stderr, "namx=%d xmin=%f bin= %f\n", nmax, xmin, bin);
  fprintf(stderr, "norm=%f\n", norm);
*/

  if( lg == 1 ) {
    to10 = 1./log(10.0);
    if(xmin <= 0.0) {
      fprintf(stderr,"min must be > 0 for -l option\n");
      exit(1);
    }
    /*    xm = log(xmin)*to10-bin/2.; */
    xm = log(xmin)*to10;
  }
  else {
    /*    xm = xmin - bin/2.;  */
    xm = xmin;
  }

  for(i=0; i<=nmax; i++){
    dn[i] =0;
    mean[i] =0;
  }




  while( scanf("%f", &x) != EOF) {
    if(lg && x <= 0.0){
    }
    else {
      xx = lg ? log(x)*to10 : x ;

      i = ( xx-xm ) / bin ; 
      /*      i = ( xx-xm ) / bin + 1; */
      /*  i>=0 is NG. since xm-bin ~ xm will give i=0 */
      if(xx >= xm && i <=nmax ) { 
      /*      if(i >= 1 && i <=nmax ) {  */
	dn[i]++;
	mean[i] += x;
      }
    }
  } 
   
  for(imin=0; dn[imin]== 0; imin++);
  for(imax=nmax; dn[imax]==0; imax--);

  sum= 0;
  for(i = imin; i<=imax; i++) {
    sum += dn[i];
  }

  /*  printf("%d\n", imin); */

  inc = lg ? pow(10.0, bin) : bin;
  /*  log(x)= log(xmin)+ imin*bin+bin/2;  
      log(x1)= log(xmin)+ imin*bin; log(x2)=log(xmin)+(imin+1)*bin;
      dx = xmin*(10**(imin+1)bin-10**(imin)*bin)
  */
  x = lg ? xmin*pow(10.0, imin*bin+bin/2) : xmin + bin*imin+bin/2;
  /*  dx = x*pow(10.0, -bin/2.0)*(inc - 1.); */
  /*  printf("%f\n", dx); */
  for(i=imin; i<=imax; i++) {
    dx =lg?  xmin*(pow(10.0, (i+1)*bin)-pow(10.0, i*bin) ): bin;
    mean[i] = dn[i] == 0 ? x : mean[i]/dn[i];
    if(norm > 0.0 && norm != 1.0 ) 
      printf("%f %10.4g %d %10.4g %f\n", x, dn[i]/dx/norm, dn[i], sum/norm, mean[i]);
    else if( norm < 0.0 ) 
      printf("%f %10.4g %d %d %f\n", x,  -dn[i]/dx/norm, dn[i], sum, mean[i]);
    else
      printf("%f %10.4g %d %d %f\n", x,  dn[i]/dx, dn[i], sum, mean[i]);
    sum -= dn[i];
    x = lg ? x*inc : x + inc;
    /*    dx  = lg ? dx*inc : bin; */
  }
}


